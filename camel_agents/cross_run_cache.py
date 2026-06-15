"""
跨运行结果缓存系统 (CrossRunCache)

避免重复计算相同结构的 SurFF、AdsorbDiff 等耗时步骤。
基于结构指纹和参数哈希实现智能缓存。
"""

import hashlib
import json
import pickle
import logging
from pathlib import Path
from typing import Any, Dict, Optional, Callable, Union
from dataclasses import dataclass
from datetime import datetime, timedelta
import sqlite3
import numpy as np

from ase import Atoms
from ase.io import read, write

logger = logging.getLogger(__name__)


class StructureFingerprinter:
    """结构指纹生成器 - 生成结构的唯一标识"""
    
    @staticmethod
    def generate_fingerprint(atoms: Atoms, tolerance: float = 0.01) -> str:
        """
        生成结构的指纹
        
        基于：
        - 化学式
        - 原子类型和相对位置（容差内）
        - 晶格参数（归一化）
        
        Returns:
            64字符的十六进制指纹
        """
        # 排序原子（按元素和位置）
        symbols = np.array(atoms.get_chemical_symbols())
        positions = atoms.get_positions()
        
        # 按元素和z坐标排序
        sort_indices = np.lexsort((positions[:, 2], positions[:, 1], positions[:, 0], symbols))
        sorted_symbols = symbols[sort_indices]
        sorted_positions = positions[sort_indices]
        
        # 量化位置（容差）
        quantized_positions = np.round(sorted_positions / tolerance).astype(int)
        
        # 构建指纹数据
        fingerprint_data = {
            "formula": atoms.get_chemical_formula(),
            "symbols": "".join(sorted_symbols),
            "positions_hash": hashlib.sha256(quantized_positions.tobytes()).hexdigest()[:32],
            "cell_hash": hashlib.sha256(np.array(atoms.cell).tobytes()).hexdigest()[:16] if atoms.cell else "no_cell",
            "pbc": tuple(bool(x) for x in atoms.pbc) if hasattr(atoms, 'pbc') else (True, True, True)
        }
        
        # 计算最终指纹
        fp_string = json.dumps(fingerprint_data, sort_keys=True)
        return hashlib.sha256(fp_string.encode()).hexdigest()
    
    @staticmethod
    def generate_params_fingerprint(params: Dict) -> str:
        """生成参数字典的指纹"""
        # 过滤掉不可序列化的参数
        serializable = {}
        for k, v in params.items():
            if isinstance(v, (int, float, str, bool, type(None))):
                serializable[k] = v
            elif isinstance(v, (list, tuple)):
                serializable[k] = list(v)
            elif isinstance(v, dict):
                serializable[k] = {kk: vv for kk, vv in v.items() 
                                 if isinstance(vv, (int, float, str, bool, type(None)))}
        
        param_string = json.dumps(serializable, sort_keys=True)
        return hashlib.md5(param_string.encode()).hexdigest()


@dataclass
class CacheEntry:
    """缓存条目"""
    cache_key: str
    operation: str  # "surff", "adsorbdiff", "neb", etc.
    structure_fingerprint: str
    params_fingerprint: str
    result_path: str
    timestamp: datetime
    access_count: int = 0
    last_accessed: Optional[datetime] = None
    
    def to_dict(self) -> Dict:
        return {
            "cache_key": self.cache_key,
            "operation": self.operation,
            "structure_fingerprint": self.structure_fingerprint,
            "params_fingerprint": self.params_fingerprint,
            "result_path": self.result_path,
            "timestamp": self.timestamp.isoformat(),
            "access_count": self.access_count,
            "last_accessed": self.last_accessed.isoformat() if self.last_accessed else None
        }
    
    @classmethod
    def from_dict(cls, data: Dict) -> "CacheEntry":
        return cls(
            cache_key=data["cache_key"],
            operation=data["operation"],
            structure_fingerprint=data["structure_fingerprint"],
            params_fingerprint=data["params_fingerprint"],
            result_path=data["result_path"],
            timestamp=datetime.fromisoformat(data["timestamp"]),
            access_count=data.get("access_count", 0),
            last_accessed=datetime.fromisoformat(data["last_accessed"]) if data.get("last_accessed") else None
        )


class CrossRunCache:
    """
    跨运行结果缓存
    
    支持 SQLite 后端，持久化存储计算结果。
    自动管理缓存大小和过期时间。
    """
    
    def __init__(
        self,
        cache_dir: str = ".catdt_cache",
        max_size_gb: float = 10.0,
        default_ttl_days: int = 30,
        enable_stats: bool = True
    ):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.max_size_gb = max_size_gb
        self.default_ttl_days = default_ttl_days
        self.enable_stats = enable_stats
        
        # 初始化数据库
        self.db_path = self.cache_dir / "cache.db"
        self._init_db()
        
        # 统计
        self.stats = {"hits": 0, "misses": 0, "saved_time": 0.0}
        
        logger.info(f"Cache initialized at {self.cache_dir}")
    
    def _init_db(self):
        """初始化 SQLite 数据库"""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS cache_entries (
                    cache_key TEXT PRIMARY KEY,
                    operation TEXT NOT NULL,
                    structure_fingerprint TEXT NOT NULL,
                    params_fingerprint TEXT NOT NULL,
                    result_path TEXT NOT NULL,
                    timestamp TEXT NOT NULL,
                    access_count INTEGER DEFAULT 0,
                    last_accessed TEXT
                )
            """)
            
            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_operation 
                ON cache_entries(operation)
            """)
            
            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_fingerprint 
                ON cache_entries(structure_fingerprint)
            """)
            
            conn.commit()
    
    def get(
        self,
        operation: str,
        structure: Atoms,
        params: Optional[Dict] = None,
        ttl_days: Optional[int] = None
    ) -> Optional[Any]:
        """
        获取缓存结果
        
        Parameters:
            operation: 操作类型 ("surff", "adsorbdiff", etc.)
            structure: 输入结构
            params: 参数字典
            ttl_days: 自定义过期时间
        
        Returns:
            缓存的结果，如果没有则返回 None
        """
        structure_fp = StructureFingerprinter.generate_fingerprint(structure)
        params_fp = StructureFingerprinter.generate_params_fingerprint(params or {})
        cache_key = f"{operation}_{structure_fp}_{params_fp}"
        
        # 查询数据库
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute(
                "SELECT result_path, timestamp, access_count FROM cache_entries WHERE cache_key = ?",
                (cache_key,)
            )
            row = cursor.fetchone()
            
            if row is None:
                self.stats["misses"] += 1
                return None
            
            result_path, timestamp_str, access_count = row
            timestamp = datetime.fromisoformat(timestamp_str)
            
            # 检查过期
            ttl = ttl_days or self.default_ttl_days
            if datetime.now() - timestamp > timedelta(days=ttl):
                logger.debug(f"Cache entry expired: {cache_key}")
                self._remove_entry(cache_key)
                self.stats["misses"] += 1
                return None
            
            # 加载结果
            try:
                result = self._load_result(result_path)
                
                # 更新访问统计
                conn.execute(
                    """UPDATE cache_entries 
                       SET access_count = ?, last_accessed = ? 
                       WHERE cache_key = ?""",
                    (access_count + 1, datetime.now().isoformat(), cache_key)
                )
                conn.commit()
                
                self.stats["hits"] += 1
                logger.debug(f"Cache hit: {operation}")
                return result
                
            except Exception as e:
                logger.warning(f"Failed to load cached result: {e}")
                self._remove_entry(cache_key)
                return None
    
    def set(
        self,
        operation: str,
        structure: Atoms,
        params: Optional[Dict],
        result: Any,
        metadata: Optional[Dict] = None
    ) -> str:
        """
        保存结果到缓存
        
        Returns:
            缓存键
        """
        structure_fp = StructureFingerprinter.generate_fingerprint(structure)
        params_fp = StructureFingerprinter.generate_params_fingerprint(params or {})
        cache_key = f"{operation}_{structure_fp}_{params_fp}"
        
        # 保存结果文件
        result_path = self._save_result(cache_key, result)
        
        # 保存到数据库
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(
                """INSERT OR REPLACE INTO cache_entries 
                   (cache_key, operation, structure_fingerprint, params_fingerprint, 
                    result_path, timestamp, access_count, last_accessed)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                (cache_key, operation, structure_fp, params_fp, result_path,
                 datetime.now().isoformat(), 0, None)
            )
            conn.commit()
        
        logger.debug(f"Cache saved: {operation}")
        
        # 检查缓存大小
        self._enforce_size_limit()
        
        return cache_key
    
    def compute_or_cache(
        self,
        operation: str,
        structure: Atoms,
        compute_fn: Callable[[], Any],
        params: Optional[Dict] = None,
        ttl_days: Optional[int] = None,
        force_recompute: bool = False
    ) -> Any:
        """
        计算或从缓存获取
        
        这是一个便捷方法，封装了 get/set 逻辑。
        
        Example:
            >>> result = cache.compute_or_cache(
            ...     "surff",
            ...     bulk_structure,
            ...     lambda: run_surff_calculation(bulk_structure),
            ...     params={"top_n": 5}
            ... )
        """
        if not force_recompute:
            cached = self.get(operation, structure, params, ttl_days)
            if cached is not None:
                logger.info(f"Using cached {operation} result")
                return cached
        
        # 执行计算
        logger.info(f"Computing {operation}...")
        start_time = datetime.now()
        result = compute_fn()
        compute_time = (datetime.now() - start_time).total_seconds()
        
        # 保存到缓存
        self.set(operation, structure, params, result)
        
        # 记录节省的时间（下次使用）
        self.stats["saved_time"] += compute_time
        
        return result
    
    def invalidate(
        self,
        operation: Optional[str] = None,
        structure: Optional[Atoms] = None
    ):
        """
        使缓存失效
        
        - operation=None, structure=None: 清空所有缓存
        - operation="surff", structure=None: 清空所有 surff 缓存
        - operation="surff", structure=bulk: 清空该结构的 surff 缓存
        """
        with sqlite3.connect(self.db_path) as conn:
            if operation is None and structure is None:
                # 清空所有
                conn.execute("DELETE FROM cache_entries")
                logger.info("Cache cleared")
            elif structure is None:
                # 清空特定操作
                conn.execute("DELETE FROM cache_entries WHERE operation = ?", (operation,))
                logger.info(f"Cache cleared for operation: {operation}")
            else:
                # 清空特定结构和操作
                structure_fp = StructureFingerprinter.generate_fingerprint(structure)
                conn.execute(
                    "DELETE FROM cache_entries WHERE operation = ? AND structure_fingerprint = ?",
                    (operation, structure_fp)
                )
                logger.info(f"Cache cleared for {operation} on structure")
            
            conn.commit()
    
    def get_stats(self) -> Dict:
        """获取缓存统计"""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("SELECT COUNT(*) FROM cache_entries")
            total_entries = cursor.fetchone()[0]
            
            cursor = conn.execute(
                "SELECT operation, COUNT(*) FROM cache_entries GROUP BY operation"
            )
            operation_counts = {row[0]: row[1] for row in cursor.fetchall()}
        
        hit_rate = 0
        total_requests = self.stats["hits"] + self.stats["misses"]
        if total_requests > 0:
            hit_rate = self.stats["hits"] / total_requests
        
        return {
            "total_entries": total_entries,
            "operation_counts": operation_counts,
            "hits": self.stats["hits"],
            "misses": self.stats["misses"],
            "hit_rate": f"{hit_rate:.1%}",
            "saved_time_hours": self.stats["saved_time"] / 3600,
            "cache_size_mb": self._get_cache_size_mb()
        }
    
    def _save_result(self, cache_key: str, result: Any) -> str:
        """保存结果到文件"""
        # 使用两级目录结构避免单个目录文件过多
        subdir = cache_key[:2]
        result_dir = self.cache_dir / "results" / subdir
        result_dir.mkdir(parents=True, exist_ok=True)
        
        result_path = result_dir / f"{cache_key}.pkl"
        
        with open(result_path, 'wb') as f:
            pickle.dump(result, f)
        
        return str(result_path)
    
    def _load_result(self, result_path: str) -> Any:
        """从文件加载结果"""
        with open(result_path, 'rb') as f:
            return pickle.load(f)
    
    def _remove_entry(self, cache_key: str):
        """删除缓存条目"""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute(
                "SELECT result_path FROM cache_entries WHERE cache_key = ?",
                (cache_key,)
            )
            row = cursor.fetchone()
            
            if row:
                result_path = Path(row[0])
                if result_path.exists():
                    result_path.unlink()
                
                conn.execute("DELETE FROM cache_entries WHERE cache_key = ?", (cache_key,))
                conn.commit()
    
    def _enforce_size_limit(self):
        """强制执行大小限制"""
        current_size = self._get_cache_size_mb()
        max_size_mb = self.max_size_gb * 1024
        
        if current_size <= max_size_mb:
            return
        
        logger.info(f"Cache size ({current_size:.1f} MB) exceeds limit ({max_size_mb:.1f} MB), cleaning...")
        
        # 删除最旧的条目
        with sqlite3.connect(self.db_path) as conn:
            while current_size > max_size_mb * 0.9:  # 清理到 90%
                cursor = conn.execute(
                    """SELECT cache_key, result_path FROM cache_entries 
                       ORDER BY COALESCE(last_accessed, timestamp) ASC LIMIT 1"""
                )
                row = cursor.fetchone()
                
                if row is None:
                    break
                
                self._remove_entry(row[0])
                current_size = self._get_cache_size_mb()
    
    def _get_cache_size_mb(self) -> float:
        """获取缓存大小（MB）"""
        total_size = 0
        for path in self.cache_dir.rglob("*"):
            if path.is_file():
                total_size += path.stat().st_size
        return total_size / (1024 * 1024)


# Forward declaration to avoid circular import
class CachedCatDTTools:
    """
    带缓存的 CatDTTools 包装器
    
    自动缓存 SurFF、AdsorbDiff 等耗时操作的结果。
    """
    
    def __init__(self, tools: Any, cache: Optional[CrossRunCache] = None):
        self.tools = tools
        self.cache = cache or CrossRunCache()
    
    def generate_surfaces(
        self,
        bulk_structure_path: str,
        top_n_surfaces: int = 3,
        run_id: str = "default_run",
        workflow_step: str = "01_surfaces",
        **kwargs
    ) -> str:
        """缓存版本的 generate_surfaces"""
        from ase.io import read
        
        bulk = read(bulk_structure_path)
        params = {"top_n": top_n_surfaces, **kwargs}
        
        def compute():
            return self.tools.generate_surfaces(
                bulk_structure_path=bulk_structure_path,
                top_n_surfaces=top_n_surfaces,
                run_id=run_id,
                workflow_step=workflow_step,
                **kwargs
            )
        
        result_path = self.cache.compute_or_cache(
            "surff", bulk, compute, params
        )
        
        # 复制结果到当前运行目录
        import shutil
        from pathlib import Path
        
        current_output = Path(self.tools.output_base_dir) / run_id / workflow_step
        current_output.mkdir(parents=True, exist_ok=True)
        
        # 读取缓存的结果并重新保存到当前位置
        result = self.tools._load_result_from_pickle(result_path)
        new_path = current_output / "surff_prediction_result.pkl"
        self.tools._save_result_to_pickle(result, new_path)
        
        return str(new_path)
    
    def predict_adsorption_sites(
        self,
        surface_path: str,
        adsorbate_smi: str,
        num_sites: int = 5,
        **kwargs
    ) -> str:
        """缓存版本的 predict_adsorption_sites"""
        from ase.io import read
        
        surface = read(surface_path)
        params = {"adsorbate": adsorbate_smi, "num_sites": num_sites, **kwargs}
        
        def compute():
            return self.tools.predict_adsorption_sites(
                surface_path=surface_path,
                adsorbate_smi=adsorbate_smi,
                num_sites=num_sites,
                **kwargs
            )
        
        result_path = self.cache.compute_or_cache(
            "adsorbdiff", surface, compute, params
        )
        
        # 复制到当前位置
        from pathlib import Path
        result = self.tools._load_result_from_pickle(result_path)
        
        run_id = kwargs.get('run_id', 'default_run')
        workflow_step = kwargs.get('workflow_step', '02_adsorption')
        current_output = Path(self.tools.output_base_dir) / run_id / workflow_step
        current_output.mkdir(parents=True, exist_ok=True)
        
        new_path = current_output / "adsorbdiff_prediction_result.pkl"
        self.tools._save_result_to_pickle(result, new_path)
        
        return str(new_path)
    
    def get_cache_stats(self) -> Dict:
        """获取缓存统计"""
        return self.cache.get_stats()


# =============================================================================
# 便捷函数
# =============================================================================

def clear_all_cache(cache_dir: str = ".catdt_cache"):
    """清空所有缓存"""
    cache = CrossRunCache(cache_dir=cache_dir)
    cache.invalidate()
    print("All cache cleared")


def show_cache_stats(cache_dir: str = ".catdt_cache"):
    """显示缓存统计"""
    cache = CrossRunCache(cache_dir=cache_dir)
    stats = cache.get_stats()
    
    print("\n" + "="*50)
    print("Cache Statistics")
    print("="*50)
    print(f"Total entries: {stats['total_entries']}")
    print(f"Cache size: {stats['cache_size_mb']:.1f} MB")
    print(f"Hit rate: {stats['hit_rate']}")
    print(f"Hits: {stats['hits']}, Misses: {stats['misses']}")
    print(f"Estimated saved time: {stats['saved_time_hours']:.1f} hours")
    print("\nBy operation:")
    for op, count in stats['operation_counts'].items():
        print(f"  {op}: {count}")


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        if sys.argv[1] == "clear":
            clear_all_cache()
        elif sys.argv[1] == "stats":
            show_cache_stats()
    else:
        print("Usage: python cross_run_cache.py [clear|stats]")

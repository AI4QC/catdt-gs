# Files

## File: camel_agents/tooling/__init__.py
````python
"""Tool modules for CatDT CAMEL integration."""
⋮----
__all__ = [
````

## File: camel_agents/tooling/care_bridge.py
````python
"""CARE bridge layer for CatDT.

Wraps CARE's blueprint expansion and template chains, converting CARE
objects into CatDT's own ``CatalyticStateRecord`` / ``ElementaryStepCandidate``
representations.  CARE is treated as a pluggable backend — if it is
unavailable or the query falls outside its domain, the bridge returns an
explicit rejection rather than crashing.
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
# ---------------------------------------------------------------------------
# Lazy CARE imports — fail gracefully
⋮----
_CARE_AVAILABLE: Optional[bool] = None
⋮----
def _check_care_available() -> bool
⋮----
import care.crn.utils.blueprint  # noqa: F401
_CARE_AVAILABLE = True
⋮----
_CARE_AVAILABLE = False
⋮----
# CARE's hard-coded element space
_CARE_SUPPORTED_ELEMENTS = frozenset({"C", "H", "O", "N"})
⋮----
# Domain report
⋮----
class CAREDomainReport(BaseModel)
⋮----
"""Whether the query falls inside CARE's supported domain."""
in_domain: bool = False
unsupported_elements: List[str] = Field(default_factory=list)
reason: str = ""
care_available: bool = False
⋮----
# Bridge class
⋮----
class CAREBridge
⋮----
"""Adapter between CARE reaction-network generation and CatDT schemas."""
⋮----
def __init__(self, num_cpu: int = 1, show_progress: bool = False)
⋮----
# ----- domain check -----
⋮----
"""Check whether the given species fall inside CARE's supported domain."""
care_available = _check_care_available()
⋮----
species_elements = set()
⋮----
# Extract elements from labels like "*CO", "H2(g)", "CH3OH"
cleaned = label.replace("*", "").replace("(g)", "").replace("(s)", "")
⋮----
found = re.findall(r"[A-Z][a-z]?", cleaned)
⋮----
# Remove metal surface elements (not part of CARE's species space)
adsorbate_elements = species_elements - {
⋮----
unsupported = sorted(adsorbate_elements - _CARE_SUPPORTED_ELEMENTS)
⋮----
# ----- blueprint generation -----
⋮----
"""Generate a CARE reaction network blueprint from seed species.

        Returns a dict with keys: ``network``, ``intermediates``, ``reactions``,
        ``stats``, or ``error`` on failure.
        """
⋮----
network = gen_blueprint(
⋮----
# ----- individual template expansion -----
⋮----
def expand_dissociation(self, chemical_space: List[str]) -> Dict[str, Any]
⋮----
"""Run CARE dissociation template on a list of species SMILES/formulas."""
⋮----
def expand_adsorption(self, intermediates_dict: Any) -> List[Dict[str, Any]]
⋮----
"""Run CARE adsorption template on intermediates."""
⋮----
reactions = gen_adsorption_reactions(
⋮----
def expand_rearrangement(self, intermediates_dict: Any) -> List[Dict[str, Any]]
⋮----
"""Run CARE rearrangement template on intermediates."""
⋮----
reactions = gen_rearrangement_reactions(
⋮----
def expand_pcet(self, intermediates_dict: Any, reactions_list: Any) -> List[Dict[str, Any]]
⋮----
"""Run CARE PCET template on intermediates + existing reactions."""
⋮----
pcet_reactions = gen_pcet_reactions(
⋮----
# ----- conversion to CatDT schemas -----
⋮----
"""Convert a CARE ReactionNetwork to CatDT state records + step candidates.

        Returns:
            (state_records, step_candidates) — both as list of dicts
            matching CatalyticStateRecord / ElementaryStepCandidate field names.
        """
state_records = []
step_candidates = []
⋮----
intermediates = self._extract_intermediates(network)
⋮----
reactions = self._extract_reactions(network)
⋮----
# ----- internal helpers -----
⋮----
@staticmethod
    def _is_intermediate(node: Any) -> bool
⋮----
def _extract_intermediates(self, network: Any) -> Dict[str, Dict[str, Any]]
⋮----
result = {}
⋮----
key = getattr(node, "code", None) or str(id(node))
⋮----
def _extract_reactions(self, network: Any) -> List[Dict[str, Any]]
⋮----
result = []
⋮----
rxn_dict = self._reaction_to_dict(node)
⋮----
@staticmethod
    def _intermediate_to_dict(inter: Any) -> Dict[str, Any]
⋮----
formula = getattr(inter, "formula", "") or ""
phase = getattr(inter, "phase", "ads") or "ads"
smiles = getattr(inter, "smiles", "") or ""
code = getattr(inter, "code", "") or ""
⋮----
# Parse element counts from formula
⋮----
elements: Dict[str, int] = {}
⋮----
elem = match.group(1)
count = int(match.group(2) or 1)
⋮----
@staticmethod
    def _reaction_to_dict(rxn: Any) -> Dict[str, Any]
⋮----
rxn_type = getattr(rxn, "r_type", "") or getattr(rxn, "reaction_type", "") or "other"
⋮----
reactants = []
products = []
reactant_key = ""
product_key = ""
bond_changes = ""
⋮----
label = getattr(r, "formula", str(r))
⋮----
reactant_key = getattr(r, "code", label)
⋮----
label = getattr(p, "formula", str(p))
⋮----
product_key = getattr(p, "code", label)
⋮----
bb = rxn.bb
bond_changes = f"break {getattr(bb, 'atom1', '?')}-{getattr(bb, 'atom2', '?')}"
⋮----
bond_changes = str(rxn.bond_changes)
⋮----
# Module-level convenience functions
⋮----
def care_check_domain(species_labels: List[str]) -> CAREDomainReport
⋮----
"""Quick domain check without instantiating the full bridge."""
⋮----
"""One-shot CRN blueprint generation."""
⋮----
__all__ = [
````

## File: camel_agents/tooling/common.py
````python
"""Shared runtime utilities for CatDT CAMEL tools."""
⋮----
CATDT_CORE_PATH = Path(__file__).resolve().parents[2] / "core"
DEPS_BASE_PATH = CATDT_CORE_PATH.parent / "deps"
⋮----
logger = logging.getLogger(__name__)
⋮----
class CatDTToolRuntimeBase
⋮----
# VSSR-MC energy model and electrochemical parameters
⋮----
# Set MP API key in environment if provided
⋮----
def _get_dt_instance(self, use_visualization: bool = False) -> Any
⋮----
surff_root = str(DEPS_BASE_PATH / "SurFF")
adsorbdiff_root = str(DEPS_BASE_PATH / "AdsorbDiff")
surface_sampling_root = str(DEPS_BASE_PATH / "surface-sampling")
fairchem_root = str(DEPS_BASE_PATH / "fairchem")
fairchem_model_path = str(DEPS_BASE_PATH / "fairchem_models/uma-s-1p1.pt")
⋮----
# Map config mc_energy_model to VSSRMCPredictor model_type
vssr_model = "CHGNetNFF"
⋮----
vssr_model = "UMA"
⋮----
shared_kwargs = dict(
⋮----
def _get_viz_manager_instance(self, current_output_dir: Path) -> Any
⋮----
def _make_pickle_safe(self, obj: Any, depth: int = 0) -> Any
⋮----
safe = {"__type__": type(obj).__name__}
⋮----
def _save_result_to_pickle(self, result: Any, path: Path) -> str
⋮----
"""Helper to save complex objects to pickle file."""
⋮----
fallback = {
⋮----
def _load_result_from_pickle(self, path: str) -> Any
⋮----
"""Helper to load complex objects from pickle file."""
⋮----
predictor = self._shared_fairchem_predictors.get(cache_key)
⋮----
predictor_work_dir = self.output_base_dir / work_subdir / cache_key
⋮----
predictor = FairchemPredictor(
````

## File: camel_agents/tooling/dashboard.py
````python
"""
Workflow Dashboard - 实时监控和可视化仪表板

提供：
1. 实时工作流状态监控
2. 步骤进度追踪
3. 资源使用监控
4. 结果可视化
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
@dataclass
class StepMetrics
⋮----
"""步骤指标"""
step_name: str
status: str  # "pending", "running", "completed", "failed", "skipped"
start_time: Optional[datetime] = None
end_time: Optional[datetime] = None
progress: float = 0.0  # 0-1
metrics: Dict[str, Any] = field(default_factory=dict)
error: Optional[str] = None
⋮----
@property
    def duration(self) -> Optional[float]
⋮----
"""运行时长（秒）"""
⋮----
def to_dict(self) -> Dict
⋮----
@dataclass
class WorkflowMetrics
⋮----
"""工作流指标"""
run_id: str
start_time: datetime
⋮----
status: str = "running"  # "running", "completed", "failed", "stopped"
steps: List[StepMetrics] = field(default_factory=list)
current_step: Optional[str] = None
⋮----
@property
    def total_duration(self) -> float
⋮----
"""总运行时长"""
end = self.end_time or datetime.now()
⋮----
@property
    def overall_progress(self) -> float
⋮----
"""整体进度"""
⋮----
class WorkflowDashboard
⋮----
"""
    工作流监控仪表板
    
    生成实时监控页面和报告
    """
⋮----
refresh_interval: int = 5  # 秒
⋮----
# 历史数据（用于趋势图）
⋮----
def start_monitoring(self)
⋮----
"""开始后台监控"""
⋮----
def stop_monitoring(self)
⋮----
"""停止监控"""
⋮----
def _background_update(self)
⋮----
"""后台更新线程"""
⋮----
"""更新步骤状态"""
# 查找或创建步骤
step = next((s for s in self.metrics.steps if s.step_name == step_name), None)
⋮----
step = StepMetrics(step_name=step_name, status=status)
⋮----
# 更新状态
old_status = step.status
⋮----
# 时间戳管理
⋮----
# 更新当前步骤
⋮----
# 更新工作流状态
⋮----
# 记录历史
⋮----
def log_energy(self, step_name: str, energy: float, **kwargs)
⋮----
"""记录能量值"""
⋮----
def generate_dashboard(self) -> Path
⋮----
"""生成监控页面"""
html = self._generate_html()
⋮----
dashboard_path = self.output_dir / "index.html"
⋮----
# 同时保存 JSON 数据
json_path = self.output_dir / "metrics.json"
⋮----
def _generate_html(self) -> str
⋮----
"""生成 HTML 页面"""
m = self.metrics
⋮----
# 步骤表格行
step_rows = []
⋮----
duration_str = f"{step.duration:.1f}s" if step.duration else "N/A"
progress_bar = self._generate_progress_bar(step.progress)
status_color = {
⋮----
metrics_str = "<br>".join(f"{k}: {v}" for k, v in step.metrics.items()) if step.metrics else "-"
⋮----
step_table = "\n".join(step_rows) if step_rows else "<tr><td colspan='6'>No steps yet</td></tr>"
⋮----
# 整体进度
overall_progress_bar = self._generate_progress_bar(m.overall_progress, width=400)
⋮----
html = f"""
⋮----
def _generate_progress_bar(self, progress: float, width: int = 200) -> str
⋮----
"""生成进度条 HTML"""
percentage = min(100, max(0, progress * 100))
⋮----
def generate_summary_report(self) -> Path
⋮----
"""生成最终摘要报告"""
⋮----
report_lines = [
⋮----
duration = f"{step.duration:.1f}s" if step.duration else "N/A"
⋮----
report_text = "\n".join(report_lines)
⋮----
report_path = self.output_dir / "summary_report.txt"
⋮----
class WorkflowMonitor
⋮----
"""
    工作流监控器 - 便捷的上下文管理器
    
    Example:
        >>> with WorkflowMonitor("my_run") as monitor:
        ...     monitor.update_step("agent1", "running")
        ...     # 执行步骤
        ...     monitor.update_step("agent1", "completed", progress=1.0)
    """
⋮----
def __enter__(self)
⋮----
def __exit__(self, exc_type, exc_val, exc_tb)
⋮----
# =============================================================================
# 便捷函数
⋮----
def create_dashboard_for_run(run_id: str, output_dir: str = "output/dashboard") -> WorkflowDashboard
⋮----
"""
    为运行创建仪表板
    
    Example:
        >>> dashboard = create_dashboard_for_run("co_to_ch4_001")
        >>> dashboard.update_step("agent1_surface", "running")
        >>> # ... 执行步骤
        >>> dashboard.update_step("agent1_surface", "completed", progress=1.0)
        >>> dashboard.generate_dashboard()
    """
⋮----
def generate_workflow_summary(metrics_file: str) -> str
⋮----
"""从指标文件生成摘要"""
⋮----
data = json.load(f)
⋮----
lines = [
````

## File: camel_agents/tooling/geometry.py
````python
"""Position calculations, distance metrics, atom manipulation, PBC handling, and staging."""
⋮----
_current_file_dir = Path(__file__).parent.resolve()
_project_root = _current_file_dir.parent.parent
⋮----
class WorkflowGeometryMixin
⋮----
@staticmethod
    def _safe_position(candidate: Any, fallback: List[float], structure: Optional[Atoms] = None) -> List[float]
⋮----
position = [float(candidate[0]), float(candidate[1]), float(candidate[2])]
⋮----
cell = np.array(structure.cell)
⋮----
lengths = np.linalg.norm(cell, axis=1)
max_len = float(np.max(lengths)) if lengths.size else 0.0
⋮----
frac = np.linalg.pinv(cell.T) @ np.array(position, dtype=float)
⋮----
wrapped = cell.T @ frac
⋮----
@staticmethod
    def _resolve_atom_index(idx_value: Any, adsorbate_indices: List[int], total_atoms: int) -> int
⋮----
_ = total_atoms  # Pathway atom operations are restricted to adsorbate atoms only.
⋮----
idx = int(idx_value)
⋮----
@staticmethod
    def _minimum_distance(symbol_a: str, symbol_b: str, scale: float = 0.85) -> float
⋮----
z_a = atomic_numbers.get(symbol_a)
z_b = atomic_numbers.get(symbol_b)
⋮----
def _resolve_addition_elements(self, atom_spec: AtomAddSpec) -> List[str]
⋮----
candidates = [
⋮----
tokens: List[str] = []
⋮----
text = str(raw or "").strip()
⋮----
parsed = self._extract_species_tokens(text)
⋮----
canon = self._canonical_element(text)
⋮----
def _step_formula_target_counts(self, step: PathwayStepSpec) -> Counter
⋮----
reactant_tokens = Counter(self._extract_species_tokens(step.reactant_formula))
product_tokens = Counter(self._extract_species_tokens(step.product_formula))
⋮----
target: Counter = Counter()
⋮----
@staticmethod
    def _count_elements_on_indices(structure: Atoms, indices: List[int]) -> Counter
⋮----
counts: Counter = Counter()
⋮----
to_delete = sorted({idx for idx in delete_indices if 0 <= idx < len(structure)}, reverse=True)
⋮----
delete_set = set(to_delete)
⋮----
def _remap(old_idx: int) -> Optional[int]
⋮----
shift = sum(1 for d in to_delete if d < old_idx)
new_idx = old_idx - shift
⋮----
new_ads: List[int] = []
⋮----
mapped = _remap(int(idx))
⋮----
new_staged: List[int] = []
⋮----
ads = sorted({idx for idx in adsorbate_indices if 0 <= idx < len(structure)})
staged = sorted({idx for idx in staged_indices if 0 <= idx < len(structure)})
⋮----
removable = [idx for idx in ads if idx >= max(0, fixed_count) or idx in set(staged)]
removable_set = set(removable)
core_refs = [idx for idx in ads if idx not in removable_set]
⋮----
positions = structure.get_positions()
⋮----
anchor = np.mean(positions[core_refs], axis=0)
⋮----
anchor = np.mean(positions[ads], axis=0)
⋮----
current_counts = self._count_elements_on_indices(structure, ads)
remove_plan: List[int] = []
⋮----
target = int(target_counts.get(element, 0))
excess = int(current - target)
⋮----
candidates = [idx for idx in removable if structure[idx].symbol == element]
⋮----
add_offset = 0
⋮----
missing = int(target - current_counts.get(element, 0))
⋮----
pos = self._get_unbonded_position(structure, ads)
⋮----
new_idx = len(structure) - 1
⋮----
def _relax_close_contacts(self, structure: Atoms, movable_indices: List[int], max_iter: int = 8) -> None
⋮----
symbols = structure.get_chemical_symbols()
movable = set(movable_indices)
⋮----
changed = False
⋮----
min_dist = self._minimum_distance(symbols[i], symbols[j])
dist = float(np.linalg.norm(positions[j] - positions[i]))
⋮----
direction = positions[move_idx] - positions[anchor_idx]
norm = float(np.linalg.norm(direction))
⋮----
direction = np.array([0.0, 0.0, 1.0])
norm = 1.0
direction = direction / norm
⋮----
changed = True
⋮----
def _get_pathway_energy_helper(self) -> Optional[PathwayPredictor]
⋮----
fairchem_root = _project_root / "deps" / "fairchem"
⋮----
corrected = float(product_energy - reactant_energy)
helper = self._get_pathway_energy_helper()
⋮----
element_change = helper.calculate_element_change(reactant_label, product_label)
⋮----
def _min_pair_distance(self, structure: Atoms, indices: List[int]) -> float
⋮----
min_dist = float("inf")
⋮----
dist = float(np.linalg.norm(positions[idx_i] - positions[idx_j]))
⋮----
min_dist = dist
⋮----
@staticmethod
    def _min_distance_any_pair(structure: Atoms) -> float
⋮----
dist = float(np.linalg.norm(positions[i] - positions[j]))
⋮----
@staticmethod
    def _min_adsorbate_surface_distance(structure: Atoms, adsorbate_indices: List[int]) -> float
⋮----
valid_ads = [idx for idx in adsorbate_indices if 0 <= idx < len(structure)]
⋮----
ads_set = set(valid_ads)
surface_indices = [i for i in range(len(structure)) if i not in ads_set]
⋮----
dist = float(np.linalg.norm(positions[a_idx] - positions[s_idx]))
⋮----
frame_count = max(3, int(n_images))
reactant_pos = np.array(reactant.get_positions(), dtype=float)
product_pos = np.array(product.get_positions(), dtype=float)
frame_atoms = reactant.copy()
⋮----
min_any_pair = float("inf")
min_any_pair_frame = -1
min_ads_surface = float("inf")
min_ads_surface_frame = -1
⋮----
alpha = frame_idx / float(frame_count - 1)
interp_pos = (1.0 - alpha) * reactant_pos + alpha * product_pos
⋮----
min_any = self._min_distance_any_pair(frame_atoms)
⋮----
min_any_pair = min_any
min_any_pair_frame = frame_idx
⋮----
min_ads_surf = self._min_adsorbate_surface_distance(frame_atoms, adsorbate_indices)
⋮----
min_ads_surface = min_ads_surf
min_ads_surface_frame = frame_idx
⋮----
min_any_pair = -1.0
⋮----
min_ads_surface = -1.0
⋮----
unwrapped = np.array([positions[idx] for idx in adsorbate_indices], dtype=float)
⋮----
cell = np.array(structure.cell, dtype=float)
⋮----
inv_t = np.linalg.pinv(cell.T)
⋮----
ref_pool = [idx for idx in (reference_indices or adsorbate_indices) if idx in set(adsorbate_indices)]
⋮----
ref_pool = list(adsorbate_indices)
ref_idx = ref_pool[0]
ref_frac = inv_t @ np.array(positions[ref_idx], dtype=float)
⋮----
frac = inv_t @ np.array(positions[atom_idx], dtype=float)
delta = frac - ref_frac
⋮----
frac_aligned = ref_frac + delta
⋮----
aligned = self._unwrap_adsorbate_positions(
⋮----
"""Align only movable adsorbate atoms to core reference periodic image."""
⋮----
ordered_ads = [idx for idx in adsorbate_indices if 0 <= idx < len(structure)]
⋮----
row_map = {idx: row for row, idx in enumerate(ordered_ads)}
⋮----
row = row_map.get(idx)
⋮----
def _adsorbate_centroid_xy(self, structure: Atoms, adsorbate_indices: List[int]) -> Optional[np.ndarray]
⋮----
positions = self._unwrap_adsorbate_positions(structure, adsorbate_indices)
⋮----
def _estimate_surface_top_z(self, structure: Atoms, adsorbate_indices: List[int]) -> float
⋮----
"""Estimate top-surface z robustly from non-adsorbate dominant slab species."""
⋮----
ads_set = set(adsorbate_indices or [])
surface_candidates = [i for i in range(len(structure)) if i not in ads_set]
⋮----
dominant = Counter(symbols[i] for i in surface_candidates).most_common(1)[0][0]
dominant_indices = [i for i in surface_candidates if symbols[i] == dominant]
ref_indices = dominant_indices if dominant_indices else surface_candidates
⋮----
@staticmethod
    def _contains_gas_phase_hint(*texts: Any) -> bool
⋮----
t = str(text or "").lower()
⋮----
r_pos = reactant.get_positions()
p_pos = product.get_positions()
r_sym = reactant.get_chemical_symbols()
p_sym = product.get_chemical_symbols()
⋮----
pairs: List[Tuple[int, int]] = []
used_product: set[int] = set()
⋮----
from scipy.optimize import linear_sum_assignment  # type: ignore
⋮----
linear_sum_assignment = None
⋮----
by_elem: Dict[str, Tuple[List[int], List[int]]] = {}
⋮----
cost = np.zeros((len(r_list), len(p_list)), dtype=float)
⋮----
r_idx = r_list[ri]
p_idx = p_list[cj]
⋮----
best_idx = None
best_dist = float("inf")
⋮----
dist = float(np.linalg.norm(r_pos[r_idx] - p_pos[p_idx]))
⋮----
best_dist = dist
best_idx = p_idx
⋮----
reactant_symbols = reactant.get_chemical_symbols()
product_symbols = product.get_chemical_symbols()
⋮----
product_positions = product.get_positions()
reactant_positions = reactant.get_positions()
⋮----
assignment: Dict[int, int] = {}
by_symbol: Dict[str, Tuple[List[int], List[int]]] = {}
⋮----
cost = np.zeros((len(r_group), len(p_group)), dtype=float)
⋮----
remaining = list(p_group)
⋮----
chosen = min(
⋮----
reorder_sequence = [assignment[idx] for idx in range(len(reactant_symbols))]
⋮----
reordered_product = product[reorder_sequence]
old_to_new = {old: new for new, old in enumerate(reorder_sequence)}
new_ads = [old_to_new[idx] for idx in product_adsorbate_indices if idx in old_to_new]
new_ads = sorted(new_ads)
new_ads_set = set(new_ads)
⋮----
def _get_unbonded_position(self, structure: Atoms, adsorbate_indices: List[int]) -> List[float]
⋮----
ref = np.mean(positions[adsorbate_indices], axis=0)
surface_z = self._estimate_surface_top_z(structure, adsorbate_indices)
z_target = float(np.clip(max(ref[2] + 0.6, surface_z + 1.0), surface_z + 0.9, surface_z + 3.2))
⋮----
ref = positions[np.argmax(positions[:, 2])]
surface_z = float(np.max(positions[:, 2]))
⋮----
def _get_bonded_position(self, structure: Atoms, adsorbate_indices: List[int], element: str) -> List[float]
⋮----
anchor_idx = min(adsorbate_indices, key=lambda idx: positions[idx][2])
ref = positions[anchor_idx]
anchor_symbol = symbols[anchor_idx]
⋮----
anchor_idx = int(np.argmax(positions[:, 2]))
⋮----
bond_dist = max(1.0, self._minimum_distance(anchor_symbol, element, scale=0.9))
⋮----
"""Find a hollow/bridge adsorption site on the surface, away from adsorbate.

        Algorithm:
        1. Compute adsorbate centroid (heavy atoms: C, O, N preferred)
        2. Find top-layer surface metal atoms
        3. Group by distance from adsorbate centroid
        4. Skip nearest neighbors (< min_dist), take second-nearest group
        5. Compute geometric center of 2-3 metal atoms as hollow site
        6. z = site_atoms_z + 1.7 Å (typical adsorption height above local site)

        Returns None if no suitable site is found.
        """
⋮----
# Adsorbate centroid (heavy atoms preferred)
⋮----
heavy_ads = [i for i in adsorbate_indices if symbols[i] not in ("H",) and 0 <= i < len(structure)]
ref_ads = heavy_ads if heavy_ads else [i for i in adsorbate_indices if 0 <= i < len(structure)]
⋮----
centroid_xy = np.mean(positions[ref_ads, :2], axis=0)
⋮----
# Top-layer surface metal atoms
⋮----
surf_symbols = [symbols[i] for i in surface_indices]
dominant_elem = Counter(surf_symbols).most_common(1)[0][0]
dominant_indices = [i for i in surface_indices if symbols[i] == dominant_elem]
⋮----
# Find the dominant terrace z (most common z-level for surface atoms)
# This avoids using VSSR-MC added atoms which sit above the main terrace
dominant_z_values = [positions[i, 2] for i in dominant_indices]
⋮----
# Round to 0.1 Å bins and find the most populated z-layer
z_bins: Dict[float, List[int]] = {}
⋮----
z_bin = round(positions[idx, 2], 1)
⋮----
# Sort by population (most atoms) then by z (HIGHEST preferred = top surface terrace)
sorted_bins = sorted(z_bins.items(), key=lambda x: (-len(x[1]), -x[0]))
terrace_z = sorted_bins[0][0]
terrace_indices = sorted_bins[0][1]
⋮----
# Use terrace atoms for site finding (not VSSR-MC adatoms)
top_layer = terrace_indices
⋮----
top_layer = dominant_indices
⋮----
# Sort by xy distance from adsorbate centroid
dist_list = []
⋮----
xy_dist = float(np.linalg.norm(positions[idx, :2] - centroid_xy))
⋮----
# Find atoms in target distance range
candidates = [(idx, d) for idx, d in dist_list if min_dist_from_ads <= d <= max_dist_from_ads]
⋮----
# Widen search: 2.5-6.5 Å
candidates = [(idx, d) for idx, d in dist_list if 2.5 <= d <= 6.5]
⋮----
# Last resort: anything > 2.0 Å
candidates = [(idx, d) for idx, d in dist_list if d > 2.0]
⋮----
# Apply offset to pick DIFFERENT site groups for multiple staged atoms
# With offset=0, pick the first group; offset=1, jump ahead by group_size
group_size = min(3, len(candidates))
# For multiple staged atoms, ensure ≥2.5 Å separation by using different site groups
start = (offset * max(group_size, 2)) % len(candidates)
site_atoms = []
⋮----
site_center = np.mean([positions[idx, :2] for idx, _ in site_atoms], axis=0)
# Use the z of the actual site atoms (terrace level), not surface_top_z
site_z = float(np.mean([positions[idx, 2] for idx, _ in site_atoms]))
⋮----
z = float(site_z + 1.7 + 0.2 * offset)
⋮----
"""Enumerate ALL surface adsorption sites near a given anchor atom.

        Uses pymatgen ``AdsorbateSiteFinder`` for robust site identification on
        any surface type (metals, oxides, alloys, stepped surfaces), then filters
        by distance from the anchor atom.

        Parameters
        ----------
        structure : Atoms
            Slab+adsorbate structure.
        adsorbate_indices : list of int
            Indices of the main adsorbate atoms.
        anchor_pos : ndarray, shape (3,)
            3D position of the anchor atom (the atom the co-adsorbate will bond
            to / dissociate from in the other endpoint).
        element : str
            Element to place (determines adsorption height adjustment).
        min_dist, max_dist : float
            Distance range from anchor_pos (Å).

        Returns
        -------
        list of dict
            Each entry: {"position": [x,y,z], "dist_from_anchor": float, "site_label": str}
            Sorted by distance from anchor (closest first).
        """
subprocess_sites = self._enumerate_surface_sites_with_pymatgen_subprocess(
⋮----
@staticmethod
    def _wrap_cartesian_position(structure: Atoms, position: np.ndarray) -> np.ndarray
⋮----
@staticmethod
    def _pbc_cartesian_delta(structure: Atoms, cart_a: np.ndarray, cart_b: np.ndarray) -> np.ndarray
⋮----
delta = np.array(cart_b, dtype=float) - np.array(cart_a, dtype=float)
⋮----
frac_delta = inv_t @ delta
⋮----
def _pbc_distance(self, structure: Atoms, cart_a: np.ndarray, cart_b: np.ndarray) -> float
⋮----
slab_indices = [i for i in range(len(structure)) if i not in ads_set]
⋮----
slab_ase = structure[slab_indices].copy()
⋮----
timeout_s = max(5, int(os.getenv("CATDT_PYMATGEN_SITE_TIMEOUT_SEC", "20")))
retry_count = max(1, int(os.getenv("CATDT_PYMATGEN_SITE_RETRIES", "3")))
env = os.environ.copy()
⋮----
# Prevent child process from touching GPU — avoids CUDA fork SIGSEGV
⋮----
# Sync + release CUDA context before fork to avoid nondeterministic
# SIGSEGV in child processes (CUDA driver atexit handlers)
⋮----
slab_path = Path(tmpdir) / "slab.traj"
⋮----
payload = {
worker_path = Path(__file__).with_name("pymatgen_site_worker.py")
proc: Optional[subprocess.CompletedProcess[str]] = None
⋮----
proc = subprocess.run(
⋮----
proc = None
⋮----
stderr_tail = (proc.stderr or "").strip().splitlines()[-1] if proc.stderr else ""
⋮----
stdout = (proc.stdout or "").strip()
⋮----
parsed = json.loads(stdout)
⋮----
sites: List[Dict[str, Any]] = []
⋮----
pos = item.get("position")
⋮----
pair_distances: List[float] = []
⋮----
dxy = float(np.linalg.norm(self._pbc_cartesian_delta(structure, positions[idx_i], positions[idx_j])[:2]))
⋮----
nearest_neighbor = min(pair_distances) if pair_distances else 2.6
bridge_cutoff = nearest_neighbor * 1.35
square_diag_cutoff = nearest_neighbor * 1.70
height_map = {"H": 1.0, "O": 1.5, "C": 1.7, "N": 1.5, "S": 1.8}
site_height = float(terrace_z + height_map.get(element, max(1.2, self._minimum_distance(dominant_elem, element, scale=1.0))))
⋮----
candidates: List[Dict[str, Any]] = []
seen_xy: set[Tuple[float, float]] = set()
⋮----
def add_candidate(raw_pos: np.ndarray, label: str) -> None
⋮----
wrapped = self._wrap_cartesian_position(structure, np.array(raw_pos, dtype=float))
candidate = np.array([wrapped[0], wrapped[1], site_height], dtype=float)
dist = self._pbc_distance(structure, np.array(anchor_pos, dtype=float), candidate)
key = (round(float(candidate[0]), 2), round(float(candidate[1]), 2))
⋮----
pos_i = positions[idx_i]
⋮----
pos_j = positions[idx_j]
delta = self._pbc_cartesian_delta(structure, pos_i, pos_j)
dxy = float(np.linalg.norm(delta[:2]))
⋮----
midpoint = pos_i + 0.5 * delta
⋮----
filtered = [site for site in candidates if min_dist <= float(site["dist_from_anchor"]) <= max_dist]
⋮----
filtered = [site for site in candidates if 2.0 <= float(site["dist_from_anchor"]) <= 6.0]
⋮----
"""Find a surface adsorption site for co-adsorbate near an anchor atom.

        Parameters
        ----------
        anchor_pos : ndarray or None
            If given, sites are ranked by distance to this point (the atom the
            co-adsorbate bonds/dissociates from). If None, uses adsorbate centroid.
        """
⋮----
heavy = [i for i in adsorbate_indices
ref = heavy if heavy else [i for i in adsorbate_indices if 0 <= i < len(structure)]
⋮----
anchor_pos = np.mean(positions[ref], axis=0)
⋮----
sites = self.enumerate_surface_sites_near_anchor(
⋮----
site = sites[offset % len(sites)]
result_positions = [site["position"]]
⋮----
# Multi-atom fragments: subsequent atoms at bond distance
⋮----
prev = result_positions[ei - 1]
elem = coadsorbate_elements[ei]
bond_d = max(0.96, self._minimum_distance(
angle = ei * 1.047
dx = bond_d * 0.3 * np.cos(angle)
dy = bond_d * 0.3 * np.sin(angle)
dz = np.sqrt(max(0.1, bond_d**2 - dx**2 - dy**2))
⋮----
"""Position for an incoming (added-to-reactant) staged atom.

        Incoming atoms exist in the product but not the reactant.
        Strategy: use the OTHER endpoint (product) to find which anchor atom
        *element* bonds to, then place the staged atom above that anchor
        on THIS side (reactant) at bond distance — representing "approaching".

        This avoids interpolation collisions (staged atom stays on the same
        side as in the product) and avoids identical endpoints (staged atom
        is offset above the anchor, not at the final bonded position).
        """
⋮----
anchor_pos = self._find_bonding_anchor_on_this_side(
⋮----
"""Position for an outgoing (added-to-product) staged atom.

        Outgoing atoms exist in the reactant but not the product.
        Strategy: use the OTHER endpoint (reactant) to find which anchor atom
        *element* was bonded to, then find that anchor on THIS side (product)
        and place the staged atom displaced outward — representing "departing".
        """
⋮----
# Shift outward from adsorbate centroid (departing direction)
pos_arr = np.array(anchor_pos, dtype=float)
valid_ads = [i for i in adsorbate_indices if 0 <= i < len(structure)]
⋮----
centroid = np.mean(structure.get_positions()[valid_ads], axis=0)
direction = pos_arr - centroid
norm = np.linalg.norm(direction)
⋮----
# Extra 1.5 Å outward to distinguish from bonded position
pos_arr = pos_arr + direction * 1.5
⋮----
"""Find the anchor atom for *element* using the other endpoint, then
        return a position above the corresponding anchor on THIS side.

        1. In other_endpoint, find *element* among adsorbate atoms.
        2. Find which heavy atom it's nearest to (= bonding partner / anchor).
        3. Find the same-element anchor on this_endpoint.
        4. Place staged atom at bond distance above that anchor on this side.
        """
other_syms = other_endpoint.get_chemical_symbols()
other_pos = other_endpoint.get_positions()
⋮----
# Find *element* in other endpoint
elem_indices = [i for i in other_ads_indices
⋮----
elem_idx = elem_indices[min(offset, len(elem_indices) - 1)]
elem_pos = other_pos[elem_idx]
⋮----
# Find its nearest heavy-atom bonding partner in other endpoint
heavy_other = [i for i in other_ads_indices
⋮----
anchor_other_idx = min(heavy_other, key=lambda i: float(np.linalg.norm(other_pos[i] - elem_pos)))
anchor_sym = other_syms[anchor_other_idx]
⋮----
# Find same-element anchor on THIS side
this_syms = this_endpoint.get_chemical_symbols()
this_pos = this_endpoint.get_positions()
anchor_candidates = [i for i in this_ads_indices
⋮----
# Pick the anchor closest to surface (lowest z) — most likely the same atom
anchor_this_idx = min(anchor_candidates, key=lambda i: this_pos[i, 2])
anchor_this_pos = this_pos[anchor_this_idx]
⋮----
# Place staged atom at bond distance above this anchor
bond_dist = max(1.0, self._minimum_distance(anchor_sym, element, scale=1.0))
angle = offset * 2.094
dx = 0.3 * np.cos(angle)
dy = 0.3 * np.sin(angle)
dz = np.sqrt(max(0, bond_dist**2 - dx**2 - dy**2))
⋮----
"""Find the position of *element* among adsorbate atoms in an endpoint."""
symbols = endpoint.get_chemical_symbols()
positions = endpoint.get_positions()
matches = [i for i in ads_indices
⋮----
# Pick the offset-th match (for multiple atoms of same element)
idx = matches[min(offset, len(matches) - 1)]
⋮----
"""Place atom at bond distance from the best anchor (fallback method)."""
⋮----
heavy_ads = [i for i in adsorbate_indices
ref_ads = heavy_ads if heavy_ads else [i for i in adsorbate_indices
⋮----
anchor_idx = min(
anchor_pos = positions[anchor_idx]
anchor_sym = symbols[anchor_idx]
⋮----
angle = offset * 2.094 + angle_shift
⋮----
def _get_remote_position(self, structure: Atoms) -> List[float]
⋮----
max_z = np.max(positions[:, 2])
center_x = float(np.mean(positions[:, 0]))
center_y = float(np.mean(positions[:, 1]))
⋮----
anchor_indices = core_adsorbate_indices or adsorbate_indices
⋮----
anchor = np.mean(positions[anchor_indices], axis=0)
⋮----
radius = 1.35 + 0.15 * (offset % 3)
theta = (np.pi / 3.0) * (offset % 6)
x = float(anchor[0] + radius * np.cos(theta))
y = float(anchor[1] + radius * np.sin(theta))
z_raw = max(anchor[2] + 0.45, surface_z + 1.0 + 0.05 * offset)
z = float(np.clip(z_raw, surface_z + 0.9, surface_z + 3.2))
⋮----
fallback = self._get_unbonded_position(structure, adsorbate_indices)
⋮----
base = np.array(base_position, dtype=float)
center = np.array(center_position, dtype=float)
⋮----
direction = base - center
⋮----
direction = np.array([1.0, 0.0, 0.3], dtype=float)
⋮----
staged = base + direction * float(radial_distance)
⋮----
valid_move = [i for i in move_indices if i in set(valid_ads)]
⋮----
ref_pool = [i for i in reference_indices if i in set(valid_ads)]
⋮----
ref_pool = [i for i in valid_ads if i not in set(valid_move)]
⋮----
ref_pool = list(valid_move)
⋮----
center = np.mean(positions[ref_pool], axis=0)
surface_z = self._estimate_surface_top_z(structure, valid_ads)
⋮----
placed = [i for i in valid_ads if i not in set(valid_move)]
⋮----
base = positions[idx_move].copy()
staged = self._stage_outward_from_center(
⋮----
close_idx = -1
close_dist = float("inf")
⋮----
dist = float(np.linalg.norm(staged - positions[idx_other]))
⋮----
close_dist = dist
close_idx = idx_other
⋮----
direction = staged - positions[close_idx]
⋮----
direction = np.array([1.0, 0.0, 0.2], dtype=float)
⋮----
staged = staged + direction / norm * (min_sep - close_dist + 0.08)
⋮----
"""Return (connected_components, max_nearest_neighbor, min_pair_distance) for staged atoms."""
⋮----
adjacency: Dict[int, set[int]] = {idx: set() for idx in staged}
nearest: Dict[int, float] = {idx: float("inf") for idx in staged}
min_pair = float("inf")
⋮----
min_pair = min(min_pair, dist)
⋮----
bond_cutoff = min(
⋮----
seen: set[int] = set()
components = 0
⋮----
stack = [idx]
⋮----
current = stack.pop()
⋮----
max_nearest = max((v for v in nearest.values() if np.isfinite(v)), default=0.0)
min_pair_dist = min_pair if np.isfinite(min_pair) else -1.0
⋮----
"""Keep staged extras geometrically coherent without moving fixed core coordinates."""
⋮----
ads_set = {idx for idx in adsorbate_indices if 0 <= idx < len(structure)}
movable = sorted(
⋮----
ref_pool = [idx for idx in core_reference_indices if idx in ads_set and idx not in set(movable)]
⋮----
ref_pool = [idx for idx in ads_set if idx not in set(movable)]
⋮----
ref_pool = list(movable)
⋮----
surface_z = self._estimate_surface_top_z(structure, list(ads_set))
⋮----
local_base = np.array(
radial = local_base[:2] - center[:2]
radial_norm = float(np.linalg.norm(radial))
⋮----
radial = np.array([1.0, 0.0], dtype=float)
radial_norm = 1.0
⋮----
idx_move = movable[0]
cur = np.array(structure.positions[idx_move], dtype=float)
lateral = float(np.linalg.norm(cur[:2] - center[:2]))
height = float(cur[2] - surface_z)
needs_restage = lateral < 1.10 or lateral > 4.20 or height < 0.80 or height > 2.60
⋮----
target = np.array(
⋮----
base_xy_radius = 2.0
first = movable[0]
first_pos = np.array(
⋮----
prev_idx = first
prev_pos = first_pos.copy()
⋮----
bond_len = max(
theta = float((np.pi / 2.5) * order)
direction = np.array([np.cos(theta), np.sin(theta), 0.22], dtype=float)
direction = direction / float(np.linalg.norm(direction))
⋮----
staged = prev_pos + direction * bond_len
⋮----
prev_idx = idx_move
prev_pos = staged
⋮----
candidates: List[Tuple[float, int, float, int]] = []
⋮----
prepared = prepare_endpoint_pair_for_interpolation(
prepared_ads = sorted(
⋮----
raw_ads = sorted({i for i in (list(reactant_ads) + list(product_ads)) if 0 <= i < len(reactant)})
⋮----
def _score(item: Tuple[float, int, float, int]) -> Tuple[float, float]
⋮----
pair_score = float(min_pair if min_pair > 0 else -1.0)
surf_score = float(min_ads_surf if min_ads_surf > 0 else -1.0)
⋮----
@staticmethod
    def _canonical_element(symbol: str) -> str
⋮----
s = str(symbol or "").strip()
⋮----
def _extract_species_tokens(self, species: str) -> List[str]
⋮----
raw = str(species or "").strip()
⋮----
token_symbols: List[str] = []
primary = self._extract_primary_adsorbate_species(raw)
source = primary or raw
⋮----
# direct element
canon = self._canonical_element(source)
⋮----
# parse formula-like token
formula = source.replace("*", "").replace("(", "").replace(")", "")
⋮----
symbol = self._canonical_element(elem)
⋮----
n = int(count) if count else 1
⋮----
inferred: List[str] = []
⋮----
# First try explicit species labels from LLM output
⋮----
tokens = self._resolve_addition_elements(atom_spec)
⋮----
# Generic fallback: infer by element-count delta between product and reactant
reactant_counts = Counter(reactant.get_chemical_symbols())
product_counts = Counter(product.get_chemical_symbols())
all_elements = set(reactant_counts.keys()) | set(product_counts.keys())
⋮----
delta = product_counts.get(element, 0) - reactant_counts.get(element, 0)
⋮----
expected_tokens = self._extract_species_tokens(formula)
⋮----
expected_counts = Counter(expected_tokens)
⋮----
centroid = np.mean(positions[adsorbate_indices], axis=0)
⋮----
selected: List[int] = []
used: set[int] = set()
⋮----
candidates = [idx for idx in adsorbate_indices if symbols[idx] == elem and idx not in used]
⋮----
keep = candidates[:count]
⋮----
"""Ensure formula-required adsorbate elements exist, using staged atoms only when needed."""
⋮----
core_ads = self._select_formula_core_adsorbate_indices(structure, valid_ads, formula)
core_set = set(core_ads)
⋮----
observed_counts = Counter(symbols[idx] for idx in core_ads if 0 <= idx < len(structure))
⋮----
missing = expected_counts[elem] - observed_counts.get(elem, 0)
⋮----
same_elem_spare = [
⋮----
staged = np.array(
⋮----
_ = state, base_structure
⋮----
def _agent5_programmatic_validation(self, step_structures: List[Dict[str, Any]]) -> Dict[str, Any]
⋮----
positions = source.get_positions()
anchor_pool = [idx for idx in source_adsorbate_indices if 0 <= idx < len(source)]
⋮----
anchor = np.mean(positions[anchor_pool], axis=0)
⋮----
anchor = np.mean(positions, axis=0)
⋮----
def _rank(idx: int) -> Tuple[int, float]
⋮----
prefer_rank = 0 if idx >= max(0, int(prefer_index_at_least)) else 1
dist_rank = float(np.linalg.norm(positions[idx] - anchor))
⋮----
best_idx = sorted(candidates, key=_rank)[0]
best_pos = positions[best_idx]
````

## File: camel_agents/tooling/mechanism.py
````python
"""Mechanism search tools mixin for CatDT.

Exposes the core mechanism search capabilities as tool methods that can be
registered as CAMEL FunctionTools.
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
class MechanismToolsMixin
⋮----
"""Mixin providing mechanism search tool methods for CatDTTools."""
⋮----
_mechanism_context: Optional[Dict[str, Any]] = None
_mechanism_care_bridge: Any = None
_mechanism_care_network: Any = None
⋮----
"""Initialize the mechanism search context.

        Returns a structured context dict and a CARE domain report.
        """
⋮----
ctx = MechanismContext(
⋮----
# Check CARE domain
bridge = CAREBridge(num_cpu=1)
all_species = [initial_state, target_state] + (known_intermediates or [])
domain_report = bridge.check_domain(species_labels=all_species)
⋮----
"""Generate candidate elementary steps from all enabled backends.

        Must call ``initialize_mechanism_context`` first.
        """
⋮----
ctx = self._mechanism_context
initial = ctx["initial_state"]
elements = parse_species_elements(initial)
⋮----
# Optionally generate CARE network
care_network = None
⋮----
bridge = self._mechanism_care_bridge
domain = bridge.check_domain(
⋮----
seeds = care_seed_species or self._infer_care_seeds(ctx)
bp_result = bridge.generate_blueprint(
care_network = bp_result.get("network")
⋮----
router = CandidateGeneratorRouter(
⋮----
candidates = router.generate_candidates(
⋮----
"""Estimate free energy for a single catalytic state."""
⋮----
surface_path = None
⋮----
surface_path = self._mechanism_context.get("surface_structure_path")
⋮----
thermal_backend = ThermalStateEnergyBackend(tools=self)
router = FreeEnergyRouter(thermal_backend=thermal_backend, default_backend=backend)
⋮----
"""Run mechanism search using the selected exploration module(s).

        Args:
            exploration_mode: One of:
                - ``"agent_guided"``: LLM/user recommended pathway → fast evaluate
                - ``"systematic"``: CRN from scratch (CARE-style beam search)
                - ``"both_sequential"``: agent-guided first, then systematic
                - ``"both_parallel"``: run both, merge results
        """
⋮----
all_pathways: List[Dict[str, Any]] = []
all_stats: Dict[str, Any] = {"exploration_mode": exploration_mode}
⋮----
# --- Module 1: Agent-guided pathway exploration ---
⋮----
guided_result = self._run_agent_guided_search(ctx, energy_backend)
⋮----
# --- Module 2: Systematic CRN exploration ---
⋮----
systematic_result = self._run_systematic_search(
⋮----
result = {
# Store for workflow retrieval after agent tool call
⋮----
# ---- Module 1: Agent-guided search ----
⋮----
"""Agent-guided pathway exploration.

        - If user provided ``known_intermediates``: use them directly as a single path.
        - If only initial + target given: call LLM to recommend multiple competing paths,
          evaluate each with UMA, prune by energy.
        """
⋮----
target = ctx["target_state"]
known = ctx.get("known_intermediates", [])
surface_path = ctx.get("surface_structure_path")
fixed_site = ctx.get("fixed_adsorption_site")
surface_id = ctx.get("surface_id", "")
surface_facet = ctx.get("surface_facet", "")
constraints = ctx.get("constraints", "")
⋮----
# User gave intermediates → single path, no LLM needed
sequences = [[initial] + [k for k in known if k != initial and k != target] + [target]]
⋮----
# No intermediates → call LLM to propose multiple competing paths
⋮----
sequences = self._llm_recommend_pathways(
⋮----
# Evaluate each pathway with UMA
⋮----
router = FreeEnergyRouter(thermal_backend=thermal_backend, default_backend=energy_backend)
generator = AgentGuidedPathwayGenerator()
⋮----
total_evaluated = 0
# Cache energies: same label on same surface = same energy
energy_cache: Dict[str, Optional[float]] = {}
⋮----
steps = generator.generate_pathway_steps(sequence)
⋮----
states_data = []
⋮----
energy = energy_cache[label]
⋮----
energy = None
⋮----
est = router.estimate_state_free_energy(
energy = est.get("free_energy_eV")
⋮----
steps_data = []
⋮----
# Compute pathway total energy change for ranking
energies = [s["free_energy_eV"] for s in states_data if s["free_energy_eV"] is not None]
max_energy = max(energies) if energies else None
⋮----
# Rank by max intermediate energy (lower = better)
⋮----
"""Call LLM (as a proper CAMEL agent) to recommend multiple competing pathways.

        Uses the same CamelWorkflowAgent infrastructure as Agent1-7.
        The agent can also call mechanism tools if needed.
        """
⋮----
surface_desc = (
prompt = TaskPrompts.mechanism_recommend_pathways(
⋮----
llm = get_camel_model_backend(temperature=0.7, max_tokens=2000)
agent = CamelWorkflowAgent(
⋮----
tools=[],  # pure reasoning, no tool calls
⋮----
output = agent.run(prompt, allow_tool_calls=False)
content = str(getattr(output, "raw", "") or "")
⋮----
pathways = self._parse_pathway_json(content, initial, target)
⋮----
"""Extract pathway JSON from LLM response text."""
⋮----
# Try to find JSON array in the response
# Look for [[...], [...], ...] pattern
⋮----
data = _json.loads(match.group())
⋮----
# Validate each pathway has at least 3 steps and starts/ends correctly
valid = []
⋮----
@staticmethod
    def _fallback_generic_pathway(initial: str, target: str) -> List[List[str]]
⋮----
"""When LLM is unavailable, return a minimal generic path."""
⋮----
r_elem = parse_species_elements(initial)
t_elem = parse_species_elements(target)
⋮----
# Simple heuristic: if target has more H than initial, it's a reduction
h_diff = t_elem.get("H", 0) - r_elem.get("H", 0)
o_diff = t_elem.get("O", 0) - r_elem.get("O", 0)
⋮----
intermediates = [initial]
⋮----
# Build a minimal path by stepwise H-addition / O-removal
current = dict(r_elem)
for _ in range(20):  # safety limit
⋮----
# Remove O first (if target has fewer O)
⋮----
label = _build_fallback_label(current)
⋮----
# Then add H (if target has more H)
⋮----
# ---- Module 2: Systematic CRN search ----
⋮----
"""Systematic CRN search: RDKit enumerate → LLM chemical filter → UMA evaluate.

        Phase 1: Pure graph search (NO UMA calls) — enumerate all reachable
                 states up to max_depth using RDKit bond operations.
        Phase 2: LLM chemical plausibility filter — select top-N most
                 promising pathways based on chemical rules.
        Phase 3: UMA energy evaluation — evaluate the LLM-filtered pathways
                 and rank by adsorption energy.
        """
⋮----
# ================================================================
# Phase 1: BFS graph enumeration (no UMA) — find ALL paths to target
⋮----
heuristic_router = FreeEnergyRouter(default_backend="heuristic")
policy = PruningPolicy(max_depth=max_depth)  # only used for root init
⋮----
engine = MechanismSearchEngine(
initial_elements = parse_species_elements(initial)
⋮----
# Full graph traversal — enumerate all paths to target
raw_path_dicts = engine.enumerate_all_paths_bfs(
⋮----
n_raw = len(raw_path_dicts)
⋮----
# Convert to intermediate label sequences
raw_sequences: List[List[str]] = []
⋮----
seq = [entry["species_label"] for entry in path]
⋮----
# Phase 2: LLM chemical plausibility filter (top N)
⋮----
filtered_sequences = self._llm_filter_pathways(
⋮----
filtered_sequences = raw_sequences
⋮----
n_filtered = len(filtered_sequences)
⋮----
# Phase 3: UMA energy evaluation on filtered pathways
⋮----
uma_router = FreeEnergyRouter(
⋮----
pathways_data: List[Dict[str, Any]] = []
⋮----
est = uma_router.estimate_state_free_energy(
⋮----
# Sort by max intermediate energy
⋮----
"""Use LLM to rank pathway sequences by chemical plausibility.

        Returns the top-N most promising sequences.
        """
⋮----
paths_text = ""
⋮----
surface_desc = f"{surface_id}({surface_facet})" if surface_facet else (surface_id or "metal surface")
⋮----
prompt = f"""You are a heterogeneous catalysis expert. Evaluate these candidate reaction pathways and select the {top_n} most chemically plausible ones.
⋮----
llm = get_camel_model_backend(temperature=0.3, max_tokens=500)
⋮----
# Parse indices from response
⋮----
indices = _json.loads(match.group())
⋮----
# Convert 1-based to 0-based, filter valid
selected = []
⋮----
i = idx - 1
⋮----
"""Export retained pathways to disk and prepare Agent4/5 handoff.

        Args:
            search_result: output from ``run_mechanism_search``
            output_base_dir: base output directory
            run_id: workflow run ID
        """
⋮----
output_base_dir = "output/catdt_workflow"
⋮----
run_id = "mechanism_search"
⋮----
builder = PathwayBundleBuilder(
⋮----
retained = [p for p in search_result.get("pathways", []) if p.get("is_retained", True)]
⋮----
manifest = builder.export_multi_pathway_bundle(retained)
⋮----
# Build Agent4/5 handoff: for each retained pathway, create step-level inputs
agent45_payloads: List[Dict[str, Any]] = []
⋮----
steps = pw.get("steps", [])
states = pw.get("states", [])
⋮----
payload = {
⋮----
shortlist_result = {
⋮----
# ---- internal helpers ----
⋮----
@staticmethod
    def _infer_care_seeds(ctx: Dict[str, Any]) -> List[str]
⋮----
"""Infer CARE seed species from mechanism context."""
seeds = set()
⋮----
cleaned = label.replace("*", "").replace("(g)", "").replace("(s)", "").strip()
⋮----
def _build_fallback_label(elements: Dict[str, int]) -> str
⋮----
"""Build *CHO style label from element dict (fallback helper)."""
parts = []
⋮----
count = elements.get(elem, 0)
⋮----
count = elements[elem]
⋮----
__all__ = ["MechanismToolsMixin"]
````

## File: camel_agents/tooling/neb.py
````python
"""NEB and endpoint validation tool mixin."""
⋮----
class NEBToolsMixin
⋮----
"""Reorder product atoms to match reactant atom ordering for NEB."""
⋮----
reactant_symbols = reactant.get_chemical_symbols()
product_symbols = product.get_chemical_symbols()
⋮----
product_positions = product.get_positions()
reactant_positions = reactant.get_positions()
⋮----
from scipy.optimize import linear_sum_assignment  # type: ignore
⋮----
linear_sum_assignment = None
⋮----
assignment: Dict[int, int] = {}
by_symbol: Dict[str, Tuple[List[int], List[int]]] = {}
⋮----
cost = np.zeros((len(r_group), len(p_group)), dtype=float)
⋮----
remaining = list(p_group)
⋮----
chosen = min(
⋮----
reorder_sequence = [assignment[idx] for idx in range(len(reactant_symbols))]
⋮----
reordered_product = product[reorder_sequence]
old_to_new = {old: new for new, old in enumerate(reorder_sequence)}
new_ads = [old_to_new[idx] for idx in product_adsorbate_indices if idx in old_to_new]
new_ads = sorted(new_ads)
new_ads_set = set(new_ads)
⋮----
def _build_barrier()
⋮----
predictor = self.get_shared_fairchem_predictor(
⋮----
results: Dict[str, Any] = {}
barrier = None
⋮----
barrier = _build_barrier()
⋮----
name = step["name"]
reactant = step["reactant"].copy()
product = step["product"].copy()
reactant_ads = list(step.get("reactant_adsorbate_indices", []))
product_ads = list(step.get("product_adsorbate_indices", []))
reactant_staged = list(step.get("reactant_staged_indices", []))
product_staged = list(step.get("product_staged_indices", []))
⋮----
product_staged = []
⋮----
reactant_ads_set = set(idx for idx in reactant_ads if 0 <= idx < len(reactant))
product_ads_set = set(idx for idx in product_ads if 0 <= idx < len(product))
reactant_fixed_surface = [idx for idx in range(len(reactant)) if idx not in reactant_ads_set]
product_fixed_surface = [idx for idx in range(len(product)) if idx not in product_ads_set]
⋮----
# Pre-NEB energy screening: skip obviously bad endpoints
⋮----
_r_tmp = reactant.copy()
_p_tmp = product.copy()
⋮----
e_r = float(_r_tmp.get_potential_energy())
e_p = float(_p_tmp.get_potential_energy())
dE = abs(e_r - e_p)
⋮----
result = barrier.predict_from_structures(
⋮----
# Clear any result references that hold CUDA tensors
⋮----
alloc = torch.cuda.memory_allocated() / 1e6
resrv = torch.cuda.memory_reserved() / 1e6
⋮----
"""
        验证 NEB 初末态结构，并提供智能修复建议

        Parameters
        ----------
        initial_structure_path : str
            初态结构文件路径（VASP格式）
        final_structure_path : str
            末态结构文件路径（VASP格式）
        expected_reaction_type : str, optional
            期望的反应类型（"hydrogenation", "dehydrogenation", "dissociation", "coupling"）
        auto_fix : bool
            是否尝试自动修复问题
        run_id : str
            运行ID
        workflow_step : str
            工作流步骤标识

        Returns
        -------
        str
            验证报告JSON文件路径
        """
output_dir = self.output_base_dir / run_id / workflow_step
⋮----
# 读取结构
initial = read(initial_structure_path)
final = read(final_structure_path)
⋮----
# 检查必要的元数据
⋮----
# 执行验证
validator = EnhancedNEBValidator(auto_fix=auto_fix, fix_attempts=3)
report = validator.validate_endpoint_pair(initial, final, expected_reaction_type)
⋮----
# 如果自动修复成功，保存修复后的结构
⋮----
fixed_initial_path = output_dir / "initial_fixed.vasp"
fixed_final_path = output_dir / "final_fixed.vasp"
# 注意：这里需要 validator 返回修复后的结构，我们需要修改 validator 来支持
⋮----
# 保存验证报告
report_path = output_dir / "validation_report.json"
⋮----
# 生成可视化（如果验证失败）
⋮----
def _detect_adsorbate_indices(self, atoms: Atoms) -> List[int]
⋮----
"""
        自动检测吸附物原子索引
        假设：吸附物是表面上方 z 坐标最高的非金属原子
        """
positions = atoms.get_positions()
symbols = atoms.get_chemical_symbols()
⋮----
# 找到表面元素（数量最多的金属）
⋮----
symbol_counts = Counter(symbols)
surface_elem = max(symbol_counts, key=symbol_counts.get)
⋮----
# 表面 z 坐标（表面原子的最大 z）
surface_z = max(positions[i, 2] for i, s in enumerate(symbols) if s == surface_elem)
⋮----
# 吸附物原子：在表面上方 > 0.5 Å 的非表面原子
adsorbate_indices = []
⋮----
"""生成验证问题的可视化"""
⋮----
viz_manager = self._get_viz_manager_instance(output_dir)
⋮----
# 高亮问题的结构图
initial_viz = output_dir / "initial_with_issues.png"
final_viz = output_dir / "final_with_issues.png"
⋮----
"""
        运行 NEB 计算，带自动预验证和修复

        这是增强版的 NEB 计算，会自动：
        1. 验证初末态结构
        2. 尝试自动修复问题
        3. 使用修复后的结构运行 NEB
        4. 如果失败，提供详细的诊断报告

        Parameters
        ----------
        initial_structure_path : str
            初态结构路径
        final_structure_path : str
            末态结构路径
        expected_reaction_type : str, optional
            期望的反应类型
        n_images : int
            NEB 图像数量
        run_id : str
            运行ID
        workflow_step : str
            工作流步骤

        Returns
        -------
        str
            NEB 结果 pickle 文件路径
        """
⋮----
# 步骤 1: 验证结构
⋮----
validation_report_path = self.validate_neb_endpoints(
⋮----
# 加载验证报告
⋮----
validation_result = json.load(f)
⋮----
# 步骤 2: 根据验证结果决定下一步
status = validation_result.get('status')
⋮----
# TODO: 使用修复后的结构
⋮----
else:  # pass
⋮----
# 步骤 3: 运行 NEB（使用现有的 run_neb_for_steps 逻辑）
⋮----
# 调用现有的 NEB 方法
neb_result_path = self.run_neb_for_steps(
⋮----
pathway_pickle_path=None,  # 直接提供结构路径
⋮----
# ---------------------------------------------------------------------------
# Agent4/5 generic geometry staging and validation tools
⋮----
class Agent45WorkflowToolsMixin
⋮----
def _agent45_geometry(self, workflow: Any) -> "Agent45GeometryTools"
⋮----
def build_steps_payload(self, workflow: Any, step_structures: List[Dict[str, Any]]) -> str
⋮----
def programmatic_validation(self, workflow: Any, step_structures: List[Dict[str, Any]]) -> Dict[str, Any]
⋮----
def get_step_element_deltas(self, workflow: Any) -> List[Dict[str, Any]]
⋮----
class Agent45GeometryTools
⋮----
def __init__(self, workflow: Any)
⋮----
"""Return nearest reference position using the workflow's canonical PBC distance."""
⋮----
refs_iter = list(reference_positions)
⋮----
pos_arr = np.array(position, dtype=float)
best_idx: Optional[int] = None
best_pos: Optional[np.ndarray] = None
best_dist = float("inf")
⋮----
ref_arr = np.array(ref_pos, dtype=float)
⋮----
dist = float(self.workflow._pbc_distance(structure, pos_arr, ref_arr))
⋮----
best_dist = dist
best_idx = idx
best_pos = ref_arr
⋮----
"""Validate Agent4's resolved position; fall back to staging plan if bad.

        Checks:
        - Position must be within *max_distance_from_ads* of nearest adsorbate atom.
        - Position must not overlap with any existing atom (< *min_distance_from_any*).
        If either check fails, the staging plan position is used instead.
        """
pos_arr = np.array(resolved_pos, dtype=float)
positions = structure.get_positions()
⋮----
# Check distance to nearest adsorbate atom
valid_ads = [i for i in adsorbate_indices if 0 <= i < len(structure)]
⋮----
# Check overlap with any existing atom
⋮----
# ------------------------------------------------------------------
# Tool-callable helpers for Agent4/5 (fully generic, no reaction-specific logic)
⋮----
def get_step_element_deltas(self) -> List[Dict[str, Any]]
⋮----
"""Return per-step element deltas: what atoms must be added to which side.

        Agent4 calls this first to know exactly what staging atoms are needed.
        Pure arithmetic from reactant/product element counts — no hardcoded rules.

        Returns
        -------
        list of dict, each with:
            step_name, reactant_formula, product_formula,
            to_add_to_reactant: list of element symbols,
            to_add_to_product: list of element symbols,
            total_staged: int
        """
owner = getattr(self.workflow, "tools", self.workflow)
baseline_steps = list(getattr(owner, "_agent45_baseline_steps", []) or [])
⋮----
result: List[Dict[str, Any]] = []
⋮----
name = str(entry.get("name", ""))
reactant = entry.get("reactant")
product = entry.get("product")
⋮----
r_cnt = Counter(reactant.get_chemical_symbols())
p_cnt = Counter(product.get_chemical_symbols())
all_elems = sorted(set(list(r_cnt) + list(p_cnt)))
⋮----
to_add_r: List[str] = []
to_add_p: List[str] = []
⋮----
diff = p_cnt.get(elem, 0) - r_cnt.get(elem, 0)
⋮----
@staticmethod
    def _canonicalize_fragment_elements(elements: List[str]) -> List[str]
⋮----
heavy = [str(elem) for elem in elements if str(elem) != "H"]
hydrogens = [str(elem) for elem in elements if str(elem) == "H"]
ordered = heavy + hydrogens
⋮----
def _site_retry_owner(self) -> Any
⋮----
owner = getattr(self.workflow, "tools", None)
⋮----
def _site_retry_state(self) -> Dict[Tuple[str, str, str], int]
⋮----
owner = self._site_retry_owner()
state = getattr(owner, "_agent45_site_retry_state", None)
⋮----
state = {}
⋮----
@staticmethod
    def _site_retry_key(step_name: str, side: str, fragment_label: str) -> Tuple[str, str, str]
⋮----
state = self._site_retry_state()
raw_rank = state.get(self._site_retry_key(step_name, side, fragment_label), 1)
⋮----
rank = int(raw_rank)
⋮----
rank = 1
⋮----
candidate_sites = list(hint.get("candidate_sites", []))
⋮----
key = self._site_retry_key(
⋮----
ranked = [dict(site) for site in candidate_sites]
⋮----
anchor = np.array(anchor_pos, dtype=float)
preferred = np.array(preferred_reference_pos, dtype=float) - anchor
⋮----
preferred_norm = float(np.linalg.norm(preferred))
⋮----
preferred = preferred / preferred_norm
⋮----
def _alignment(site: Dict[str, Any]) -> float
⋮----
pos = np.array(site.get("position", []), dtype=float)
⋮----
direction = pos - anchor
⋮----
norm = float(np.linalg.norm(direction))
⋮----
"""Identify extra fragment(s) and the anchor atom mapped onto the target endpoint."""
⋮----
src_syms = source_struct.get_chemical_symbols()
src_pos = source_struct.get_positions()
tgt_syms = target_struct.get_chemical_symbols()
tgt_pos = target_struct.get_positions()
⋮----
valid_src = [i for i in source_ads if 0 <= i < len(source_struct)]
⋮----
ads_atoms = _Atoms(
cutoffs = natural_cutoffs(ads_atoms, mult=1.2)
nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
⋮----
local_to_global = {li: gi for li, gi in enumerate(valid_src)}
adj: Dict[int, set[int]] = {li: set() for li in range(len(valid_src))}
⋮----
extra_cnt = Counter(extra_elems)
best_frag: Optional[List[int]] = None
⋮----
start_sym = ads_atoms[start_li].symbol
⋮----
remaining = Counter(extra_cnt)
frag: List[int] = []
queue = [start_li]
visited = {start_li}
⋮----
cur = queue.pop(0)
cur_sym = ads_atoms[cur].symbol
⋮----
best_frag = frag
⋮----
best_frag = []
⋮----
sym = ads_atoms[li].symbol
⋮----
frag_global = [local_to_global[li] for li in best_frag]
frag_elems = self._canonicalize_fragment_elements([src_syms[gi] for gi in frag_global])
frag_set = set(best_frag)
fragment_anchor_global = next((gi for gi in frag_global if src_syms[gi] != "H"), frag_global[0] if frag_global else None)
preferred_site_reference_pos = src_pos[fragment_anchor_global] if fragment_anchor_global is not None else None
⋮----
anchor_sym = "?"
anchor_pos_on_target = None
anchor_candidates_local: set[int] = set()
⋮----
frag_centroid = np.mean(ads_atoms.get_positions()[best_frag], axis=0)
anchor_local = min(
anchor_global = local_to_global[anchor_local]
anchor_sym = src_syms[anchor_global]
anchor_src_pos = src_pos[anchor_global]
⋮----
cands = [i for i in target_ads if 0 <= i < len(target_struct) and tgt_syms[i] == anchor_sym]
⋮----
anchor_tgt = min(cands, key=lambda i: float(np.linalg.norm(tgt_pos[i] - anchor_src_pos)))
anchor_pos_on_target = tgt_pos[anchor_tgt]
⋮----
valid_tgt = [i for i in target_ads if 0 <= i < len(target_struct)]
⋮----
anchor_pos_on_target = np.mean(tgt_pos[valid_tgt], axis=0)
⋮----
def _build_step_staging_hints(self, step_entry: Dict[str, Any]) -> List[Dict[str, Any]]
⋮----
reactant = step_entry.get("reactant")
product = step_entry.get("product")
⋮----
reactant_ads = list(step_entry.get("reactant_adsorbate_indices", []))
product_ads = list(step_entry.get("product_adsorbate_indices", []))
⋮----
to_add_to_reactant: List[str] = []
to_add_to_product: List[str] = []
⋮----
hints: List[Dict[str, Any]] = []
⋮----
frag_elems = self._canonicalize_fragment_elements(list(frag_info.get("fragment", [])))
⋮----
anchor_pos_3d = frag_info.get("anchor_pos")
candidate_sites: List[Dict[str, Any]] = []
⋮----
sites = self.workflow.enumerate_surface_sites_near_anchor(
ranked_sites = self._rank_candidate_sites_by_direction(
⋮----
position = [round(float(v), 2) for v in site["position"]]
⋮----
def _ensure_step_staging_hints(self, step_entry: Dict[str, Any]) -> List[Dict[str, Any]]
⋮----
hints = step_entry.get("staging_hints")
⋮----
hints = self._build_step_staging_hints(step_entry)
⋮----
fragment_elements = list(hint.get("fragment_elements", []))
⋮----
valid_sites: List[Dict[str, Any]] = []
valid_positions: List[np.ndarray] = []
⋮----
site_pos = np.array(site.get("position", []), dtype=float)
⋮----
rank = int(site.get("rank"))
⋮----
site_pos = list(site.get("position", []))
⋮----
"""Suggest positions for ALL staged atoms across ALL steps in one call.

        Parameters
        ----------
        requests : list of dict
            Each dict: ``{"step_name": "...", "element": "H", "side": "reactant"}``.

        Returns
        -------
        dict with "suggestions": list of per-request results, each with
        position, anchor_atom, bond_distance, reason.
        """
⋮----
step_lookup = {str(e.get("name", "")): e for e in baseline_steps}
⋮----
suggestions: List[Dict[str, Any]] = []
⋮----
step_name = str(req.get("step_name", ""))
element = str(req.get("element", ""))
side = str(req.get("side", "reactant"))
⋮----
step_entry = step_lookup.get(step_name)
⋮----
atoms = step_entry.get(side)
⋮----
ads_indices = list(step_entry.get(f"{side}_adsorbate_indices", []))
⋮----
hint = self._find_matching_staging_hint(step_entry, side, element)
⋮----
allowed_rank = self._get_allowed_site_rank_for_hint(step_entry, hint)
pos = list(candidate_sites[allowed_rank - 1]["position"])
⋮----
# Get the OTHER endpoint for cross-reference (avoids interpolation collisions)
other_side = "product" if side == "reactant" else "reactant"
other_atoms = step_entry.get(other_side)
other_ads = list(step_entry.get(f"{other_side}_adsorbate_indices", []))
other_ep = other_atoms if isinstance(other_atoms, Atoms) else None
⋮----
# Find anchor atom in the OTHER endpoint, then find co-adsorption site
anchor_pos = None
⋮----
other_syms = other_ep.get_chemical_symbols()
other_positions = other_ep.get_positions()
elem_in_other = [i for i in other_ads
⋮----
ei = elem_in_other[0]
heavy = [i for i in other_ads
⋮----
heavy = [i for i in other_ads if 0 <= i < len(other_ep) and i != ei]
⋮----
anchor_other = min(heavy, key=lambda j: float(
anchor_sym = other_syms[anchor_other]
this_syms = atoms.get_chemical_symbols()
this_pos_arr = atoms.get_positions()
cands = [i for i in ads_indices
⋮----
anchor_pos = this_pos_arr[min(cands, key=lambda i: this_pos_arr[i, 2])]
⋮----
site = self.workflow.find_nearby_coadsorption_site(
⋮----
pos = site["positions"][0]
anchor_info = f"co-adsorption {site.get('site_type', '?')} site"
⋮----
pos = self.workflow._get_incoming_staging_position(
⋮----
pos = self.workflow._get_outgoing_staging_position(
anchor_info = f"fallback anchor from {other_side} endpoint"
⋮----
"""Validate proposed atom positions for ALL steps in one call.

        Parameters
        ----------
        steps : list of dict
            Each dict: ``{"step_name": "...", "atoms": [{"species": "H", "side": "reactant", "position": [x,y,z]}]}``.

        Returns
        -------
        dict with: ok (bool), per_step (list of per-step results), suggestion (str)
        """
⋮----
per_step: List[Dict[str, Any]] = []
all_ok = True
⋮----
step_name = str(step_req.get("step_name", ""))
proposed_atoms = list(step_req.get("atoms", []))
⋮----
all_ok = False
⋮----
step_ok = True
atom_results: List[Dict[str, Any]] = []
⋮----
element = str(atom_spec.get("species", atom_spec.get("element", "?")))
side = str(atom_spec.get("side", "reactant"))
position = atom_spec.get("position", [])
⋮----
step_ok = False
⋮----
pos = np.array([float(position[0]), float(position[1]), float(position[2])])
⋮----
atoms: Optional[Atoms] = step_entry.get(side)
⋮----
# Distance to nearest adsorbate
valid_ads = [i for i in ads_indices if 0 <= i < len(atoms)]
min_ads_dist = float("inf")
nearest_ads = "?"
⋮----
idx = valid_ads[int(idx_rel)] if idx_rel is not None else valid_ads[0]
⋮----
p = positions[idx]
nearest_ads = f"{symbols[idx]} at [{p[0]:.2f},{p[1]:.2f},{p[2]:.2f}]"
⋮----
# Overlap with any atom
min_any_dist = float("inf")
nearest_any = "?"
⋮----
idx = 0
⋮----
nearest_any = f"{symbols[idx]} at [{p[0]:.2f},{p[1]:.2f},{p[2]:.2f}]"
⋮----
# Surface z
surf_indices = [i for i in range(len(atoms)) if i not in set(ads_indices)]
surf_top_z = float(np.max(positions[surf_indices, 2])) if surf_indices else 0.0
⋮----
atom_ok = True
issues: List[str] = []
⋮----
atom_ok = False
⋮----
allowed_site = candidate_sites[allowed_rank - 1] if candidate_sites else None
⋮----
nearest_rank = self._candidate_site_rank(candidate_sites, nearest_site)
⋮----
site_pos = nearest_site["position"]
⋮----
allowed_pos = allowed_site["position"]
⋮----
# --- Interpolation collision check ---
# Build temporary structures with staged atoms appended, then check
# if linear interpolation between reactant and product causes collisions.
interp_issue = ""
⋮----
reactant_base = step_entry.get("reactant")
product_base = step_entry.get("product")
⋮----
tmp_r = reactant_base.copy()
tmp_p = product_base.copy()
r_ads = list(step_entry.get("reactant_adsorbate_indices", []))
p_ads = list(step_entry.get("product_adsorbate_indices", []))
⋮----
elem = str(a.get("species", a.get("element", "")))
side = str(a.get("side", ""))
apos = a.get("position", [])
⋮----
apos_f = [float(apos[0]), float(apos[1]), float(apos[2])]
⋮----
# Only check if both endpoints have same element count
⋮----
interp_mins = self.workflow._best_interpolation_mapping_metrics(
finite = [float(v) for v in interp_mins if np.isfinite(v) and v > 0]
⋮----
min_interp = min(finite)
⋮----
interp_issue = (
⋮----
interp_issue = f"Interpolation check failed: {exc}"
⋮----
suggestion = "All positions OK." if all_ok else "Some positions need adjustment — see per_step details."
⋮----
def _build_gas_phase_product_from_current(self, current: Atoms) -> Atoms
⋮----
"""Create a synthetic desorbed product for gas-phase endpoints."""
product = current.copy()
adsorbate_indices = list(product.info.get("adsorbate_indices", []))
⋮----
positions = product.get_positions()
surface_top_z = self.workflow._estimate_surface_top_z(product, adsorbate_indices)
adsorbate_com_z = float(np.mean(positions[adsorbate_indices, 2]))
cell = np.array(product.cell, dtype=float)
c_len = float(np.linalg.norm(cell[2])) if cell.shape == (3, 3) else 0.0
⋮----
target_z = max(
⋮----
# Keep gas-phase endpoint inside the simulation cell vacuum region.
upper_bound = c_len - 1.2
target_z = min(target_z, upper_bound)
⋮----
target_z = min(upper_bound, surface_top_z + 3.2)
dz = float(target_z - adsorbate_com_z)
⋮----
fixed_entry = self.workflow._get_preplaced_entry(preplaced_products, product_label)
⋮----
fixed_product = fixed_entry["structure"].copy()
product_source = str(fixed_entry.get("source", "tool_preplacement"))
product_path = str(fixed_entry.get("path", ""))
surface_indices = list(fixed_entry.get("surface_indices", []))
adsorbate_indices = list(fixed_entry.get("adsorbate_indices", []))
⋮----
fallback_product = self._build_gas_phase_product_from_current(current)
surface_indices = list(fallback_product.info.get("surface_indices", []))
adsorbate_indices = list(fallback_product.info.get("adsorbate_indices", []))
⋮----
surface_indices = [i for i in range(len(fallback_product)) if i not in set(adsorbate_indices)]
⋮----
"""Sort structure by element (atomic number) using pymatgen for consistent NEB ordering.

        After sorting, atoms are grouped by element (H, C, O, Cu, …) with a
        deterministic sub-order (fractional z, y, x).  Adsorbate / staged roles
        are tracked via ``site.properties`` so they survive the sort.
        """
⋮----
adaptor = AseAtomsAdaptor()
struct = adaptor.get_structure(atoms)
⋮----
# Tag each site with its original index so we can rebuild ads/staged lists
ads_set = set(ads_indices)
staged_set = set(staged_indices)
⋮----
site.properties["_role"] = 2   # staged (also adsorbate)
⋮----
site.properties["_role"] = 1   # adsorbate
⋮----
site.properties["_role"] = 0   # surface
⋮----
struct.sort()  # sorts by species (atomic number), then frac coords
⋮----
sorted_atoms = adaptor.get_atoms(struct)
⋮----
new_ads: List[int] = []
new_staged: List[int] = []
⋮----
role = site.properties.get("_role", 0)
⋮----
ads_set_new = set(new_ads)
⋮----
"""Sort both NEB endpoints by element (pymatgen) for consistent ordering."""
⋮----
def build_steps_payload(self, step_structures: List[Dict[str, Any]]) -> str
⋮----
"""Build a compact per-step summary with only adsorbate coordinates.

        Redundancy reduction:
        - Step i product = Step i+1 reactant, so only product adsorbate is shown per step
          (step 1 also includes reactant adsorbate as the starting point)
        - Internal metadata (signatures, fixed atom counts) omitted
        - Compact xyz format instead of POSCAR
        """
blocks: List[str] = []
⋮----
step_name = str(step.get("name", f"step_{idx}"))
reactant = step.get("reactant")
product = step.get("product")
r_formula = step.get("reactant_formula", "")
p_formula = step.get("product_formula", "")
⋮----
# Compute element delta (what Agent4 needs to add/remove)
delta_str = "none"
⋮----
parts = []
⋮----
delta_str = " ".join(parts) if parts else "none"
⋮----
header = f"=== STEP {idx}: {step_name} ({r_formula} -> {p_formula}), element_delta={delta_str} ==="
⋮----
# Show BOTH reactant and product adsorbate for every step.
# Agent4 needs both endpoints to decide where to place new atoms.
⋮----
r_ads = list(step.get("reactant_adsorbate_indices", []))
⋮----
p_ads = list(step.get("product_adsorbate_indices", []))
⋮----
hints: List[str] = []
⋮----
side = str(hint.get("side", "reactant"))
frag_label = str(hint.get("fragment_label", ""))
anchor_sym = str(hint.get("anchor_symbol", "?"))
⋮----
allowed_rank = self._get_allowed_site_rank_for_hint(step, hint)
⋮----
p = site["position"]
⋮----
def _compact_adsorbate_xyz(self, atoms: Atoms, ads_indices: List[int], label: str) -> str
⋮----
"""Compact adsorbate coordinate listing: one line per atom."""
lines = [f"  [{label}_adsorbate] indices={ads_indices}"]
⋮----
sym = atoms[i].symbol
⋮----
def build_agent5_structures_payload(self, state: WorkflowState) -> str
⋮----
def _energy_gate_owner(self) -> Any
⋮----
def _get_uma_predictor_for_energy_gate(self, output_dir: str = "") -> Any
⋮----
owner = self._energy_gate_owner()
shared_getter = getattr(owner, "get_shared_fairchem_predictor", None)
⋮----
predictor = getattr(owner, "_agent45_uma_predictor", None)
⋮----
predictor_work_dir = Path(output_dir).resolve() if output_dir else Path(getattr(owner, "output_base_dir", ".")).resolve()
predictor_work_dir = predictor_work_dir / "_agent45_energy_gate"
⋮----
predictor = FairchemPredictor(
⋮----
"""Return ONLY staged atom indices for relaxation.

        Baseline adsorbate core atoms are NOT relaxed — their positions are
        determined by tool preplacement and must remain immutable so that
        NEB endpoints reflect the intended reaction path.
        Only staged atoms (added by Agent4 for element balancing) are relaxed.
        """
staged_set = {idx for idx in staged_indices if isinstance(idx, int) and 0 <= idx < total_atoms}
⋮----
relaxed = atoms.copy()
⋮----
movable_set = set(movable_indices)
fixed = [idx for idx in range(len(relaxed)) if idx not in movable_set]
⋮----
optimizer = FIRE(relaxed, logfile=None)
ran = True
⋮----
result = predictor.predict_energy(
forces = np.array(result.forces, dtype=float)
force_norm = np.linalg.norm(forces, axis=1) if forces.size else np.zeros((len(atoms),), dtype=float)
⋮----
valid_ads = [idx for idx in adsorbate_indices if 0 <= idx < len(atoms)]
⋮----
ads_forces = force_norm[valid_ads]
max_force_ads = float(np.max(ads_forces))
mean_force_ads = float(np.mean(ads_forces))
⋮----
max_force_ads = float(np.max(force_norm)) if len(force_norm) else 0.0
mean_force_ads = float(np.mean(force_norm)) if len(force_norm) else 0.0
⋮----
valid_focus = [
⋮----
focus_forces = force_norm[valid_focus]
max_force_focus = float(np.max(focus_forces))
mean_force_focus = float(np.mean(focus_forces))
⋮----
max_force_focus = max_force_ads
mean_force_focus = mean_force_ads
⋮----
"""Agent4/5 endpoint UMA energy gate: evaluate energies and forces.

        No relaxation — structures are used as-is from Agent4.
        Only computes energies, forces, and overlap diagnostics.
        """
⋮----
predictor = self._get_uma_predictor_for_energy_gate(output_dir=output_dir)
⋮----
fatal_issues: List[str] = []
warning_issues: List[str] = []
summaries: List[str] = []
step_reports: List[Dict[str, Any]] = []
⋮----
step_name = str(step.get("name", f"step_{step_idx}"))
step_report: Dict[str, Any] = {"step_index": step_idx, "step_name": step_name, "endpoints": {}}
⋮----
atoms = step[endpoint]
adsorbate_indices = list(step.get(f"{endpoint}_adsorbate_indices", []))
⋮----
result = self._evaluate_endpoint_with_uma(
⋮----
# Overlap check
min_pair = result["min_pair_distance"]
⋮----
# Force check (warning only, no fatal — NEB handles high forces)
⋮----
status = "FAIL" if fatal_issues else "PASS"
report = {
⋮----
gate_path = Path(output_dir) / "agent45_energy_gate.json"
⋮----
intermediates = self.workflow._get_working_intermediate_sequence(state)
⋮----
baseline_steps: List[Dict[str, Any]] = []
current = base_structure.copy()
⋮----
reactant_label = intermediates[idx]
product_label = intermediates[idx + 1]
⋮----
reactant = current.copy()
product = fixed_product.copy()
⋮----
reactant_ads = list(reactant.info.get("adsorbate_indices", []))
product_ads = list(product.info.get("adsorbate_indices", []))
⋮----
step_name = f"{reactant_label}_to_{product_label}"
entry = {
⋮----
current = fixed_product.copy()
⋮----
"""Build one step endpoints while preserving immutable core coordinates.

        Agent4 must provide explicit atoms_to_add with positions. No hardcoded fallback.
        """
reactant_neb = current_physical.copy()
product_neb = fixed_product.copy()
⋮----
reactant_ads = list(reactant_neb.info.get("adsorbate_indices", []))
product_ads = list(product_neb.info.get("adsorbate_indices", []))
⋮----
reactant_fixed_len = len(reactant_neb)
product_fixed_len = len(product_neb)
reactant_fixed_snapshot = np.array(reactant_neb.positions[:reactant_fixed_len], dtype=float)
product_fixed_snapshot = np.array(product_neb.positions[:product_fixed_len], dtype=float)
⋮----
reactant_staged_indices: List[int] = []
product_staged_indices: List[int] = []
⋮----
has_llm_additions = bool(step.atoms_to_add)
base_reactant_counts = Counter(reactant_neb.get_chemical_symbols())
base_product_counts = Counter(product_neb.get_chemical_symbols())
pending_staging: List[Tuple[str, int]] = []
⋮----
delta = int(base_product_counts.get(element, 0) - base_reactant_counts.get(element, 0))
⋮----
llm_positions: Dict[str, Dict[str, List[List[float]]]] = {"reactant": {}, "product": {}}
⋮----
def _append_llm_position(side: str, element: str, value: Any) -> None
⋮----
pos = [float(value[0]), float(value[1]), float(value[2])]
⋮----
elements = self.workflow._resolve_addition_elements(atom_spec)
⋮----
base_candidate = atom_spec.position
reactant_candidate = atom_spec.reactant_position if atom_spec.reactant_position is not None else base_candidate
product_candidate = atom_spec.product_position if atom_spec.product_position is not None else base_candidate
⋮----
incoming_offsets: Counter = Counter()
outgoing_offsets: Counter = Counter()
⋮----
def _pop_position(source: Dict[str, Dict[str, List[List[float]]]], side: str, element: str) -> Optional[List[float]]
⋮----
slots = source.get(side, {}).get(element, [])
⋮----
reactant_candidate = (
⋮----
reactant_pos = self.workflow._safe_position(
⋮----
product_candidate = (
⋮----
product_pos = self.workflow._safe_position(
⋮----
ignored_additions = sum(
⋮----
remove_indices = self.workflow._resolve_removal_indices(step, reactant_ads)
⋮----
illegal = [idx for idx in remove_indices if idx < reactant_fixed_len]
⋮----
removed_records: List[Dict[str, Any]] = []
⋮----
fallback = self.workflow._get_unbonded_position(product_neb, product_ads)
⋮----
position = self.workflow._safe_position(
⋮----
# Agent4 positions are used as-is. No alignment or geometric relaxation.
# Agent5 will validate and provide feedback if positions are problematic.
⋮----
reactant_surface = [i for i in range(len(reactant_neb)) if i not in set(reactant_ads)]
product_surface = [i for i in range(len(product_neb)) if i not in set(product_ads)]
⋮----
step_structures: List[Dict[str, Any]] = []
baseline_lookup = {
⋮----
built = self.apply_step_ops_on_fixed_product(
⋮----
reactant_neb = built["reactant"]
product_neb = built["product"]
⋮----
# Fix adsorbate drift: if product adsorbate centroid xy is far from
# reactant, shift product adsorbate xy only (preserve z to avoid
# embedding atoms into the surface for species with different
# adsorption heights, e.g. CH3 on surface vs CH4 in gas phase).
r_ads_idx = list(built["reactant_adsorbate_indices"])
p_ads_idx = list(built["product_adsorbate_indices"])
⋮----
r_centroid = np.mean([reactant_neb.positions[i] for i in r_ads_idx if i < len(reactant_neb)], axis=0)
p_centroid = np.mean([product_neb.positions[i] for i in p_ads_idx if i < len(product_neb)], axis=0)
xy_drift = float(np.linalg.norm(r_centroid[:2] - p_centroid[:2]))
⋮----
shift_xy = r_centroid[:2] - p_centroid[:2]
⋮----
# Reorder product atoms to match reactant ordering so that
# interpolation and NEB see consistent atom indices.
product_ads_built = list(built["product_adsorbate_indices"])
⋮----
product_neb = reordered_prod
⋮----
built["product_staged_indices"] = []  # indices invalidated by reorder
⋮----
baseline_entry = baseline_lookup.get(step.step_name, {})
baseline_hints = list(baseline_entry.get("staging_hints", [])) if isinstance(baseline_entry, dict) else []
copied_hints = [
⋮----
"""Strict no-op postprocess: verify immutable signatures, do not alter coordinates."""
processed: List[Dict[str, Any]] = []
⋮----
reactant: Atoms = entry["reactant"]
product: Atoms = entry["product"]
⋮----
r_fixed = int(entry.get("reactant_fixed_atom_count", len(reactant)))
p_fixed = int(entry.get("product_fixed_atom_count", len(product)))
⋮----
r_ref = str(entry.get("reactant_fixed_reference_signature", "")).strip()
p_ref = str(entry.get("product_fixed_reference_signature", "")).strip()
⋮----
r_now = self.workflow._fixed_coordinate_signature(reactant, r_fixed)
p_now = self.workflow._fixed_coordinate_signature(product, p_fixed)
⋮----
processed_entry = dict(entry)
⋮----
"""Generate actionable suggestion for the closest overlapping atom pair."""
⋮----
n = len(atoms)
⋮----
min_d = float("inf")
⋮----
d = float(np.linalg.norm(positions[i] - positions[j]))
⋮----
min_d = d
⋮----
pi = positions[best_i]
pj = positions[best_j]
dz_needed = max(0.5, threshold - min_d + 0.3)
higher = best_j if pj[2] >= pi[2] else best_i
elem_h = symbols[higher]
ph = positions[higher]
⋮----
def programmatic_validation(self, step_structures: List[Dict[str, Any]]) -> Dict[str, Any]
⋮----
"""Generic geometry gate for Agent5."""
⋮----
first = step_structures[0]
⋮----
def _core_signature_order_invariant(atoms: Atoms, fixed_count: int, decimals: int = 8) -> tuple
⋮----
n = max(0, min(int(fixed_count), len(atoms)))
⋮----
rows = [
⋮----
previous_product: Atoms | None = None
previous_product_fixed = 0
⋮----
reactant: Atoms = step["reactant"]
product: Atoms = step["product"]
⋮----
reactant_fixed = int(step.get("reactant_fixed_atom_count", len(reactant)))
product_fixed = int(step.get("product_fixed_atom_count", len(product)))
⋮----
reactant_ref_sig = str(step.get("reactant_fixed_reference_signature", "")).strip()
product_ref_sig = str(step.get("product_fixed_reference_signature", "")).strip()
reactant_sig = self.workflow._fixed_coordinate_signature(reactant, reactant_fixed)
product_sig = self.workflow._fixed_coordinate_signature(product, product_fixed)
⋮----
reactant_label = self.workflow._canonical_species_label(str(step.get("reactant_formula", "")))
product_label = self.workflow._canonical_species_label(str(step.get("product_formula", "")))
add_count = int(step.get("atoms_to_add_count", 0) or 0)
remove_count = int(step.get("atoms_to_remove_count", 0) or 0)
raw_reactant_formula = str(step.get("reactant_formula", "")).strip()
raw_product_formula = str(step.get("product_formula", "")).strip()
phase_change_only = (
⋮----
endpoints_identical = (
⋮----
# Gas-phase adsorption/desorption: identical coordinates are expected
# at this stage; the actual gas-phase displacement is applied later
# by _build_gas_phase_product_from_current.
⋮----
endpoint_min = min(
⋮----
# H-X bonds are shorter (O-H ~0.97, C-H ~1.09, N-H ~1.01)
# so use a lower threshold when H is involved in the closest pair
h_involved = False
⋮----
dists = struct.get_all_distances(mic=True)
⋮----
sym = struct.get_chemical_symbols()
⋮----
h_involved = True
⋮----
fatal_threshold = 0.50 if h_involved else 0.65
warning_threshold = 0.70 if h_involved else 0.80
⋮----
overlap_suggestion = self._overlap_suggestion(reactant, product, step_name, fatal_threshold)
⋮----
overlap_suggestion = self._overlap_suggestion(reactant, product, step_name, warning_threshold)
⋮----
# Check staged atoms distance from adsorbate
# Fatal if >5.0 Å (disconnected), warning if >4.0 Å
# Staged atoms at 2-4 Å are expected (approaching/departing species)
⋮----
staged_ep = list(step.get(staged_key, []))
⋮----
core_ads = [i for i in ads_ep if i not in set(staged_ep) and 0 <= i < len(atoms_ep)]
⋮----
positions_ep = atoms_ep.get_positions()
symbols_ep = atoms_ep.get_chemical_symbols()
⋮----
sp = positions_ep[si]
⋮----
nearest_idx = core_ads[int(idx_rel)]
nearest_sym = symbols_ep[nearest_idx]
⋮----
nearest_pos = positions_ep[nearest_idx]
suggest_pos = nearest_pos.copy()
⋮----
suggest_str = (
⋮----
endpoint_atoms = reactant if side == "reactant" else product
staged_indices = list(step.get(f"{side}_staged_indices", []))
⋮----
anchor_element = str(fragment_elements[0])
matching_staged = [
⋮----
endpoint_positions = endpoint_atoms.get_positions()
best_idx = min(
⋮----
sp = endpoint_positions[best_idx]
⋮----
allowed_site = candidate_sites[allowed_rank - 1]
⋮----
finite_interp = [float(v) for v in interp_mins if np.isfinite(v) and v > 0]
⋮----
min_interp = min(finite_interp)
⋮----
compare_n = min(previous_product_fixed, reactant_fixed, len(previous_product), len(reactant))
⋮----
previous_sig = _core_signature_order_invariant(previous_product, compare_n)
current_sig = _core_signature_order_invariant(reactant, compare_n)
⋮----
previous_product = product.copy()
previous_product_fixed = product_fixed
⋮----
feedback_parts: List[str] = []
⋮----
fatal_text = "\n".join(fatal_issues)
retry_site_lines: List[str] = []
⋮----
step_name = str(step.get("name", ""))
⋮----
endpoint_atoms = step.get(side)
⋮----
current_rank = self._candidate_site_rank(candidate_sites, nearest_site)
⋮----
next_rank = min(current_rank + 1, len(candidate_sites))
⋮----
next_site = candidate_sites[next_rank - 1]
line = (
⋮----
site3 = candidate_sites[next_rank]
````

## File: camel_agents/tooling/persistence.py
````python
"""
Workflow persistence: checkpoint system + reporting tools.

Provides:
1. Atomic checkpoint save/load for workflow state
2. Resumable workflow mixin for checkpoint-based restart
3. Step structure export and HTML report generation
"""
⋮----
# ---------------------------------------------------------------------------
# Data classes
⋮----
@dataclass
class StepArtifact
⋮----
"""步骤产物"""
name: str
file_path: str
file_type: str  # "pickle", "vasp", "json", etc.
description: str = ""
⋮----
def load(self) -> Any
⋮----
"""加载产物"""
path = Path(self.file_path)
⋮----
@dataclass
class WorkflowCheckpoint
⋮----
"""工作流检查点"""
run_id: str
step_name: str  # 已完成的步骤名称
step_number: int  # 步骤序号
timestamp: datetime
state_data: Dict[str, Any]  # 序列化的状态数据
artifacts: List[StepArtifact]  # 步骤产物
metadata: Dict[str, Any] = field(default_factory=dict)
⋮----
def to_dict(self) -> Dict
⋮----
@classmethod
    def from_dict(cls, data: Dict) -> "WorkflowCheckpoint"
⋮----
"""从字典恢复"""
artifacts = [
⋮----
# Checkpoint manager
⋮----
class CheckpointManager
⋮----
"""
    检查点管理器

    管理检查点的保存、加载和清理
    """
⋮----
def __init__(self, checkpoint_dir: str = ".catdt_checkpoints")
⋮----
state: Any,  # WorkflowState 对象
⋮----
"""
        原子性保存检查点

        先写入临时文件，再重命名，确保不会出现损坏的检查点
        """
timestamp = datetime.now()
⋮----
# 序列化状态
state_data = self._serialize_state(state)
⋮----
# 保存产物
artifact_list = []
run_checkpoint_dir = self.checkpoint_dir / run_id
⋮----
artifact_path = self._save_artifact(
⋮----
# 创建检查点
checkpoint = WorkflowCheckpoint(
⋮----
# 原子性保存
checkpoint_file = run_checkpoint_dir / f"checkpoint_{step_number:02d}_{step_name}.pkl"
temp_file = checkpoint_file.with_suffix(".tmp")
⋮----
# 原子性重命名
⋮----
# 同时保存 JSON 版本便于查看
json_file = checkpoint_file.with_suffix(".json")
⋮----
# 清理临时文件
⋮----
def load_latest_checkpoint(self, run_id: str) -> Optional[WorkflowCheckpoint]
⋮----
"""加载最新的检查点"""
run_dir = self.checkpoint_dir / run_id
⋮----
# 找到最新的检查点文件
checkpoint_files = sorted(
⋮----
latest_file = checkpoint_files[0]
⋮----
checkpoint = pickle.load(f)
⋮----
def load_checkpoint_by_step(self, run_id: str, step_name: str) -> Optional[WorkflowCheckpoint]
⋮----
"""按步骤名称加载检查点"""
⋮----
matches = list(run_dir.glob(f"checkpoint_*_{step_name}.pkl"))
⋮----
def list_checkpoints(self, run_id: str) -> List[Dict]
⋮----
"""列出所有检查点"""
⋮----
checkpoints = []
⋮----
data = json.load(f)
⋮----
def clean_old_checkpoints(self, run_id: str, keep_last: int = 3)
⋮----
"""清理旧检查点，只保留最近的几个"""
⋮----
# 同时删除对应的 JSON 文件
json_file = old_file.with_suffix(".json")
⋮----
def _serialize_state(self, state: Any) -> Dict[str, Any]
⋮----
"""序列化状态对象"""
⋮----
def _serialize_value(self, value: Any) -> Any
⋮----
"""序列化单个值"""
⋮----
# ASE Atoms 特殊处理
⋮----
# 尝试 pickle
⋮----
def _is_serializable(self, value: Any) -> bool
⋮----
"""检查值是否可序列化"""
⋮----
"""保存产物"""
artifact_dir = run_dir / f"step_{step_number:02d}_artifacts"
⋮----
# 根据数据类型选择保存方式
⋮----
file_path = artifact_dir / f"{name}.vasp"
⋮----
file_path = artifact_dir / f"{name}.json"
⋮----
# 默认 pickle
file_path = artifact_dir / f"{name}.pkl"
⋮----
# Resumable workflow mixin
⋮----
class ResumableWorkflowMixin
⋮----
"""
    可恢复工作流的 Mixin 类

    添加到工作流类中以支持断点续传
    """
⋮----
def __init__(self, *args, enable_checkpointing: bool = True, **kwargs)
⋮----
def run_with_resume(self, config: Any, start_from: Optional[str] = None) -> Any
⋮----
"""
        支持断点续传的运行

        Parameters
        ----------
        config : Any
            工作流配置
        start_from : str, optional
            从指定步骤开始（如果提供，则忽略检查点）
        """
run_id = getattr(config, 'run_id', 'default_run')
⋮----
# 确定起始步骤
⋮----
start_step = start_from
state = None
⋮----
checkpoint = self.checkpoint_manager.load_latest_checkpoint(run_id)
⋮----
start_step = self._get_next_step(checkpoint.step_name)
state = self._restore_state_from_checkpoint(checkpoint)
⋮----
start_step = self._step_names[0] if self._step_names else "start"
⋮----
# 执行剩余步骤
remaining_steps = self._get_steps_from(start_step)
⋮----
# 执行步骤
state = self._execute_step(step_name, state, config)
⋮----
# 保存检查点
⋮----
def _save_step_checkpoint(self, run_id: str, step_name: str, state: Any)
⋮----
"""保存步骤检查点"""
⋮----
artifacts = self._collect_step_artifacts(step_name, state)
metadata = {
⋮----
"step_duration": None,  # 可以添加计时
⋮----
def _get_next_step(self, completed_step: str) -> str
⋮----
"""获取下一步"""
⋮----
idx = self._step_names.index(completed_step)
⋮----
def _get_steps_from(self, start_step: str) -> List[str]
⋮----
"""获取从某步骤开始的所有步骤"""
⋮----
idx = self._step_names.index(start_step)
⋮----
def _restore_state_from_checkpoint(self, checkpoint: WorkflowCheckpoint) -> Any
⋮----
"""从检查点恢复状态 - 子类应重写此方法"""
⋮----
def _execute_step(self, step_name: str, state: Any, config: Any) -> Any
⋮----
"""执行单个步骤 - 子类应重写此方法"""
⋮----
def _collect_step_artifacts(self, step_name: str, state: Any) -> Dict[str, Any]
⋮----
"""收集步骤产物 - 子类可重写此方法"""
⋮----
def with_checkpoint(step_name: str, step_number: int)
⋮----
"""
    装饰器：自动为函数保存检查点

    Example:
        >>> @with_checkpoint("surff", 1)
        ... def run_surff(state, config):
        ...     # 执行步骤
        ...     return new_state
    """
def decorator(func: Callable)
⋮----
def wrapper(self, *args, **kwargs)
⋮----
# 执行函数
result = func(self, *args, **kwargs)
⋮----
run_id = getattr(kwargs.get('config'), 'run_id', 'default_run')
⋮----
# Agent4/5 Memento-style case memory (zero-parameter, memory-only tuning)
⋮----
class Agent45MementoToolsMixin
⋮----
"""
    Memento-style memory utilities for Agent4/5.

    Design goals:
    1) Zero-parameter adaptation (no model weight updates).
    2) Persist successful/failed trajectories as retrievable cases.
    3) Provide positive/negative in-context examples to improve next iterations.
    """
⋮----
def _agent45_memento_casebank_path(self) -> Path
⋮----
path = self.output_base_dir / "agent45_memento_casebank.jsonl"
⋮----
def _agent45_knowledge_bank_path(self) -> Path
⋮----
override = str(os.getenv("CATDT_AGENT45_KNOWLEDGE_BANK_PATH", "")).strip()
⋮----
path = self.output_base_dir / "agent45_knowledge_bank.jsonl"
⋮----
def _agent45_skill_bank_path(self) -> Path
⋮----
override = str(os.getenv("CATDT_AGENT45_SKILL_BANK_PATH", "")).strip()
⋮----
path = self.output_base_dir / "agent45_skill_bank.jsonl"
⋮----
@staticmethod
    def _agent45_memento_env_flag(name: str, default: str = "0") -> bool
⋮----
value = str(os.getenv(name, default)).strip().lower()
⋮----
def _load_agent45_memento_cases(self) -> List[Dict[str, Any]]
⋮----
path = self._agent45_memento_casebank_path()
⋮----
def _load_agent45_jsonl_items(self, path: Path) -> List[Dict[str, Any]]
⋮----
items: List[Dict[str, Any]] = []
⋮----
raw = line.strip()
⋮----
obj = json.loads(raw)
⋮----
def _load_agent45_knowledge_items(self) -> List[Dict[str, Any]]
⋮----
def _load_agent45_skill_items(self) -> List[Dict[str, Any]]
⋮----
"""
        Reuse Memento's pair extraction if available; otherwise use compatible fallback.
        """
np_memory_path = Path(__file__).resolve().parents[2] / "deps" / "Memento" / "memory" / "np_memory.py"
⋮----
spec = importlib.util.spec_from_file_location("memento_np_memory", str(np_memory_path))
⋮----
module = importlib.util.module_from_spec(spec)
⋮----
extract_pairs = getattr(module, "extract_pairs", None)
⋮----
pairs: List[Tuple[str, Any, int]] = []
⋮----
@staticmethod
    def _tokenize_case_text(text: str) -> List[str]
⋮----
cleaned = str(text or "").strip().lower()
⋮----
@staticmethod
    def _extract_transition_tokens(text: str) -> List[str]
⋮----
raw = str(text or "")
transitions = re.findall(r"(\*?[A-Za-z0-9\(\)\+\-]+)\s*(?:->|→)\s*(\*?[A-Za-z0-9\(\)\+\-]+)", raw)
⋮----
@staticmethod
    def _agent45_memento_mode() -> str
⋮----
mode = str(os.getenv("CATDT_AGENT45_MEMENTO_MODE", "parametric")).strip().lower()
⋮----
def _agent45_parametric_model_name(self) -> str
⋮----
def _agent45_parametric_model_path(self) -> Path
⋮----
override = str(os.getenv("CATDT_AGENT45_PARAMETRIC_MODEL_PATH", "")).strip()
⋮----
@staticmethod
    def _agent45_parametric_device() -> str
⋮----
@staticmethod
    def _agent45_case_query_text(case: Dict[str, Any]) -> str
⋮----
text = str(case.get("question", "")).strip()
⋮----
reaction_type = str(case.get("reaction_type", "")).strip()
transitions = list(case.get("transition_signature", []) or [])
intermediates = list(case.get("intermediates", []) or [])
⋮----
@staticmethod
    def _agent45_case_plan_text(case: Dict[str, Any]) -> str
⋮----
raw = case.get("plan")
⋮----
payload = {
⋮----
def _agent45_case_to_icl_text(self, case: Dict[str, Any]) -> str
⋮----
case_text = self._agent45_case_query_text(case)
plan_text = self._agent45_case_plan_text(case)
parts = [f"[CASE]\n{case_text}\n[PLAN]\n{plan_text}"]
⋮----
# Add surface context for parametric retriever training
surf = case.get("surface_snapshot", {})
⋮----
top2 = surf.get("top2_layers", [])
⋮----
elem_counts = dict(Counter(a.get("elem", "?") for a in top2)) if top2 else {}
⋮----
# Add per-step geometry (compact) for retriever training
geo_parts = []
⋮----
name = snap.get("name", "")
⋮----
ep_data = snap.get(ep, {})
⋮----
staged = ep_data.get("staged", [])
⋮----
staged_str = " ".join(
⋮----
def _agent45_parametric_cache(self) -> Dict[str, Any]
⋮----
cache = getattr(self, "_agent45_parametric_retriever_cache", None)
⋮----
cache = {}
⋮----
model_path = Path(model_path or self._agent45_parametric_model_path()).expanduser()
model_name = str(model_name or self._agent45_parametric_model_name()).strip()
device = str(device or self._agent45_parametric_device()).strip().lower()
⋮----
mtime = float(model_path.stat().st_mtime)
cache = self._agent45_parametric_cache()
cache_key = (
⋮----
module_path = Path(__file__).resolve().parents[2] / "deps" / "Memento" / "memory" / "train_memory_retriever.py"
⋮----
spec = importlib.util.spec_from_file_location("memento_train_memory_retriever", str(module_path))
⋮----
model_cls = getattr(module, "MemoryRetrieverClassifier", None)
⋮----
tokenizer = AutoTokenizer.from_pretrained(pretrained_model_name_or_path=model_name, use_fast=True)
backbone = AutoModel.from_pretrained(pretrained_model_name_or_path=model_name, use_safetensors=True)
model = model_cls(backbone).to(device)
state_dict = torch.load(str(model_path), map_location=device)
⋮----
cache = self._load_agent45_parametric_retriever(
tokenizer = cache["tokenizer"]
model = cache["model"]
target_device = str(cache["device"])
⋮----
case_texts = [self._agent45_case_to_icl_text(case) for case in cases]
natural_texts = [str(query_text)] * len(case_texts)
⋮----
batch_size = max(1, int(os.getenv("CATDT_AGENT45_PARAMETRIC_BATCH_SIZE", "32")))
max_len = max(64, int(os.getenv("CATDT_AGENT45_PARAMETRIC_MAX_LEN", "256")))
⋮----
probs: List[float] = []
⋮----
end = min(start + batch_size, len(case_texts))
t1 = tokenizer(
t2 = tokenizer(
ids1 = t1["input_ids"].to(target_device)
mask1 = t1["attention_mask"].to(target_device)
ids2 = t2["input_ids"].to(target_device)
mask2 = t2["attention_mask"].to(target_device)
logits = model(ids1, mask1, ids2, mask2)
batch_probs = torch.softmax(logits, dim=1)[:, 1].detach().cpu().tolist()
⋮----
@staticmethod
    def _extract_reaction_type_from_query(query_text: str) -> str
⋮----
match = re.search(r"reaction_type=(\w+)", str(query_text or ""))
⋮----
def _score_agent45_case(self, query_text: str, case: Dict[str, Any]) -> float
⋮----
query_tokens = set(self._tokenize_case_text(query_text))
case_tokens = set(
⋮----
token_score = 0.0
⋮----
token_score = len(query_tokens & case_tokens) / max(len(query_tokens | case_tokens), 1)
⋮----
query_transitions = set(self._extract_transition_tokens(query_text))
case_transitions = set(
trans_score = 0.0
⋮----
trans_score = len(query_transitions & case_transitions) / max(len(query_transitions), 1)
⋮----
reward = float(case.get("reward", 0.0) or 0.0)
reward_bonus = 0.05 if reward > 0 else 0.0
⋮----
# Reaction type matching: boost same-type cases, penalize cross-type
query_rxn = self._extract_reaction_type_from_query(query_text)
case_rxn = str(case.get("reaction_type", "")).strip().upper()
rxn_bonus = 0.0
⋮----
rxn_bonus = 0.15 if query_rxn == case_rxn else -0.10
⋮----
score = 0.50 * token_score + 0.30 * trans_score + reward_bonus + rxn_bonus
⋮----
def _score_agent45_knowledge_or_skill_item(self, query_text: str, item: Dict[str, Any]) -> float
⋮----
item_tokens = set(
⋮----
token_score = len(query_tokens & item_tokens) / max(len(query_tokens | item_tokens), 1)
⋮----
item_transitions = set(
⋮----
trans_score = len(query_transitions & item_transitions) / max(len(query_transitions), 1)
⋮----
confidence = float(item.get("confidence", item.get("score", 0.0)) or 0.0)
confidence_bonus = 0.05 * max(0.0, min(confidence, 1.0))
⋮----
# Reaction type matching: boost same-type items, penalize cross-type
⋮----
item_rxn = str(item.get("reaction_type", "")).strip().upper()
⋮----
rxn_bonus = 0.12 if query_rxn == item_rxn else -0.08
⋮----
score = 0.55 * token_score + 0.25 * trans_score + confidence_bonus + rxn_bonus
⋮----
@staticmethod
    def _agent45_item_short_text(item: Dict[str, Any], keys: List[str], max_len: int = 220) -> str
⋮----
value = str(item.get(key, "") or "").strip().replace("\n", " ")
⋮----
lines: List[str] = []
⋮----
title = self._agent45_item_short_text(item, ["title", "name", "id"], max_len=80) or f"{bank_name}_item_{idx}"
evidence = self._agent45_item_short_text(
transitions = list(item.get("transition_signature", []) or [])
⋮----
transitions = list(item.get("transitions", []) or [])
⋮----
items = self._load_agent45_knowledge_items()
⋮----
scored: List[Dict[str, Any]] = []
⋮----
score = self._score_agent45_knowledge_or_skill_item(query_text=query_text, item=item)
⋮----
obj = dict(item)
⋮----
selected = scored[: max(1, int(top_k))]
⋮----
items = self._load_agent45_skill_items()
⋮----
cases = self._load_agent45_memento_cases()
⋮----
_ = self._memento_extract_pairs(cases, key_field="question", value_field="plan")
⋮----
score = self._score_agent45_case(query_text=query_text, case=case)
⋮----
label = str(case.get("case_label", "")).strip().lower()
⋮----
label = "positive" if float(case.get("reward", 0.0) or 0.0) > 0 else "negative"
⋮----
item = dict(case)
⋮----
selected = scored[: max(int(top_k), 0)]
⋮----
positive = [c for c in selected if str(c.get("case_label", "")).lower() == "positive"]
negative = [c for c in selected if str(c.get("case_label", "")).lower() != "positive"]
⋮----
prompt_block = self.build_agent45_memento_prompt_block(
⋮----
model_path = self._agent45_parametric_model_path().expanduser()
⋮----
probs = self._score_agent45_cases_parametric(query_text=query_text, cases=cases)
⋮----
score = float(prob)
⋮----
"""
        Retrieve Memento-like cases for Agent4/5 prompt conditioning.
        """
⋮----
mode = self._agent45_memento_mode()
⋮----
def _short(text: Any, max_len: int = 220) -> str
⋮----
s = str(text or "").strip().replace("\n", " ")
⋮----
transitions = ", ".join(list(case.get("transition_signature", []) or [])[:4])
feedback = _short(case.get("feedback", ""), 160)
⋮----
# Surface info (cell + top-2-layer atom count)
⋮----
n_top2 = len(surf.get("top2_layers", []))
⋮----
# Include adsorbate geometry from successful cases
⋮----
core = ep_data.get("core", [])
⋮----
parts = [f"{a['elem']}{a['xyz']}" for a in core]
⋮----
d = sa.get("dist_to_nearest_ads", "?")
⋮----
issues = list(case.get("issues", []) or [])
issue_preview = "; ".join(_short(x, 120) for x in issues[:2])
⋮----
# Surface info
⋮----
# Show adsorbate geometry from failed cases
⋮----
"""
        Append one memory case in Memento-compatible JSONL format.
        """
⋮----
status = str(validation_status or "").upper()
⋮----
reward = 1.0 if status == "PASS" else 0.0
⋮----
transition_signature = [str(x).strip().lower() for x in (transition_signature or []) if str(x).strip()]
query_text = (
⋮----
case_label = "positive" if float(reward) > 0 else "negative"
case = {
⋮----
# Surface top-2-layer atoms + cell (shared across all steps)
⋮----
# Per-step adsorbate geometry: core + staged atom coords
⋮----
# Memento-compatible fields
⋮----
@staticmethod
    def _agent45_transition_overlap(case_a: Dict[str, Any], case_b: Dict[str, Any]) -> float
⋮----
a = {str(x).strip().lower() for x in list(case_a.get("transition_signature", []) or []) if str(x).strip()}
b = {str(x).strip().lower() for x in list(case_b.get("transition_signature", []) or []) if str(x).strip()}
⋮----
rng = random.Random(int(seed))
normalized: List[Dict[str, Any]] = []
⋮----
label = str(item.get("case_label", "")).strip().lower()
⋮----
label = "positive" if float(item.get("reward", 0.0) or 0.0) > 0 else "negative"
⋮----
positive_cases = [c for c in normalized if str(c.get("case_label", "")) == "positive"]
⋮----
pairs: List[Dict[str, Any]] = []
⋮----
query_text = self._agent45_case_query_text(anchor)
⋮----
pos_pool: List[Dict[str, Any]] = []
neg_pool: List[Dict[str, Any]] = []
⋮----
overlap = self._agent45_transition_overlap(anchor, cand)
is_positive_case = str(cand.get("case_label", "")) == "positive"
⋮----
fallback_pos = [c for c in positive_cases if c is anchor]
⋮----
fallback_pos = [rng.choice(positive_cases)]
pos_pool = fallback_pos
⋮----
pos_selected = pos_pool[: max(1, int(max_pos_per_query))]
neg_selected = neg_pool[: max(1, int(max_neg_per_query))]
⋮----
subset = list(cases[-max(1, int(max_cases)):])
pairs = self._build_agent45_parametric_training_pairs(
⋮----
n_pos = sum(int(x.get("label", 0)) for x in pairs)
n_neg = len(pairs) - n_pos
⋮----
retriever_dir = self._agent45_parametric_model_path().expanduser().parent
⋮----
history_path = retriever_dir / "train_history.jsonl"
⋮----
split_idx = int(round(len(pairs) * (1.0 - float(val_ratio))))
split_idx = max(1, min(split_idx, len(pairs) - 1))
train_pairs = pairs[:split_idx]
val_pairs = pairs[split_idx:]
⋮----
class _PairDataset(Dataset)
⋮----
def __init__(self, rows: List[Dict[str, Any]])
⋮----
def __len__(self) -> int
⋮----
def __getitem__(self, idx: int) -> Dict[str, Any]
⋮----
loader_generator = torch.Generator()
⋮----
def _collate(rows: List[Dict[str, Any]]) -> Tuple[Any, Any, Any, Any, Any]
⋮----
case_texts = [str(r["icl"]) for r in rows]
query_texts = [str(r["query"]) for r in rows]
labels = torch.tensor([int(r["label"]) for r in rows], dtype=torch.long)
t1 = tokenizer(case_texts, padding=True, truncation=True, max_length=max_len, return_tensors="pt")
t2 = tokenizer(query_texts, padding=True, truncation=True, max_length=max_len, return_tensors="pt")
⋮----
train_loader = DataLoader(
val_loader = DataLoader(
⋮----
optimizer = torch.optim.AdamW(model.parameters(), lr=float(learning_rate))
class_weight = torch.tensor(
criterion = nn.CrossEntropyLoss(weight=class_weight)
⋮----
best_val_acc = -1.0
best_state: Optional[Dict[str, Any]] = None
train_loss_curve: List[float] = []
val_acc_curve: List[float] = []
val_loss_curve: List[float] = []
⋮----
epoch_losses: List[float] = []
⋮----
ids1 = ids1.to(device)
mask1 = mask1.to(device)
ids2 = ids2.to(device)
mask2 = mask2.to(device)
labels = labels.to(device)
⋮----
loss = criterion(logits, labels)
⋮----
correct = 0
total = 0
val_loss_items: List[float] = []
⋮----
vloss = criterion(logits, labels)
⋮----
pred = torch.argmax(logits, dim=1)
⋮----
val_acc = float(correct / max(total, 1))
val_loss = float(sum(val_loss_items) / max(len(val_loss_items), 1))
⋮----
best_val_acc = val_acc
best_state = {k: v.detach().cpu() for k, v in model.state_dict().items()}
⋮----
meta = {
⋮----
# Invalidate cached retriever to ensure new checkpoint is loaded on next read.
⋮----
def get_agent45_parametric_retriever_status(self) -> Dict[str, Any]
⋮----
meta_path = model_path.parent / "train_meta.json"
payload: Dict[str, Any] = {
⋮----
# Reporting tools mixin
⋮----
class ReportingToolsMixin
⋮----
"""Export per-step reactant/product structures for audit before/after LLM edits."""
out_dir = Path(output_dir)
⋮----
visualizer = None
⋮----
visualizer = CatalystSurfaceVisualizer(renderer="tachyon", quality="low", auto_expand=False)
⋮----
manifest: Dict[str, Any] = {"prefix": prefix, "steps": []}
⋮----
step_name = str(step.get("name", f"step_{idx}"))
safe_name = ''.join(ch if (ch.isalnum() or ch in ('_', '-')) else '_' for ch in step_name)
⋮----
reactant = step.get("reactant")
product = step.get("product")
step_entry: Dict[str, Any] = {
⋮----
reactant_path = out_dir / f"{prefix}_step{idx:02d}_{safe_name}_reactant.vasp"
⋮----
reactant_png = out_dir / f"{prefix}_step{idx:02d}_{safe_name}_reactant.png"
⋮----
product_path = out_dir / f"{prefix}_step{idx:02d}_{safe_name}_product.vasp"
⋮----
product_png = out_dir / f"{prefix}_step{idx:02d}_{safe_name}_product.png"
⋮----
manifest_path = out_dir / f"{prefix}_manifest.json"
⋮----
"""
        Generates a comprehensive HTML summary report of the entire CatDT workflow.

        Args:
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "06_report").
            surff_result_path (str, optional): Path to SurFFPredictionResult pickle.
            adsorbdiff_result_path (str, optional): Path to AdsorbDiffPredictionOutput pickle.
            vssr_mc_result_path (str, optional): Path to VSSRMCResult pickle.
            pathway_result_path (str, optional): Path to CompletePathwayResult pickle.
            kmc_result_path (str, optional): Path to CatMAPResult pickle.
            llm_review_results (Dict[str, Any], optional): Dictionary of LLM review outputs for various steps.

        Returns:
            str: Path to the generated HTML summary file.
        """
output_dir = self.output_base_dir / run_id / workflow_step
⋮----
viz_manager = self._get_viz_manager_instance(output_dir)
all_visualizations = {
⋮----
# Collect existing visualization paths from previous steps
surff_viz_dir = self.output_base_dir / run_id / "01_surfaces" / "visualizations"
⋮----
ads_viz_dir = self.output_base_dir / run_id / "02_adsorption" / "visualizations"
⋮----
mc_viz_dir = self.output_base_dir / run_id / "03_reconstruction" / "visualizations"
⋮----
pathway_viz_dir = self.output_base_dir / run_id / "04_pathway_analysis" / "visualizations"
⋮----
kmc_viz_dir = self.output_base_dir / run_id / "05_kmc_simulation" / "visualizations"
⋮----
# Load actual result objects if paths are provided
full_results = {}
⋮----
html_path = viz_manager.generate_summary_html(
⋮----
def check_checkpoint_exists(self, run_id: str) -> Optional[str]
⋮----
"""
        检查是否存在可恢复的检查点

        Returns
        -------
        str or None
            最新的检查点步骤名称，如果没有则返回 None
        """
checkpoint_manager = CheckpointManager()
checkpoint = checkpoint_manager.load_latest_checkpoint(run_id)
⋮----
def list_available_checkpoints(self, run_id: str) -> List[Dict]
⋮----
"""
        列出所有可用的检查点

        Returns
        -------
        List[Dict]
            检查点列表，每个包含 step, timestamp 等信息
        """
````

## File: camel_agents/tooling/pymatgen_site_worker.py
````python
"""Isolated pymatgen adsorption-site enumeration worker.

This module is executed in a subprocess so that pymatgen/ruamel import crashes
cannot take down the main CatDT workflow process.
"""
⋮----
def main() -> int
⋮----
payload = json.loads(sys.argv[1])
slab = read(payload["slab_path"])
anchor = payload["anchor_pos"]
min_dist = float(payload["min_dist"])
max_dist = float(payload["max_dist"])
⋮----
pmg_slab = AseAtomsAdaptor.get_structure(slab)
asf = AdsorbateSiteFinder(pmg_slab)
raw_sites = asf.find_adsorption_sites(symm_reduce=0)
lattice = Lattice(pmg_slab.lattice.matrix)
⋮----
def pbc_distance(cart_a, cart_b)
⋮----
frac_a = lattice.get_fractional_coords(cart_a)
frac_b = lattice.get_fractional_coords(cart_b)
⋮----
sites = []
seen_xy = set()
⋮----
dist = pbc_distance(pos3d, anchor)
⋮----
key = (round(float(pos3d[0]), 1), round(float(pos3d[1]), 1))
````

## File: camel_agents/tooling/simulation.py
````python
"""Surface/reconstruction/pathway/KMC simulation tool mixin."""
⋮----
VisualizedGasSolidDigitalTwin = type(None)
⋮----
class SimulationToolsMixin
⋮----
"""
        Generates the most exposed surfaces from a bulk crystal structure using SurFF.

        Args:
            bulk_structure_path (str): Absolute path to the bulk crystal structure file (e.g., POSCAR).
            top_n_surfaces (int): Number of top most exposed surfaces to generate and analyze.
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "01_surfaces").

        Returns:
            str: Path to a pickle file containing the SurFFPredictionResult object,
                 and the output directory where individual slab files are saved.
        """
output_dir = self.output_base_dir / run_id / workflow_step
⋮----
dt = self._get_dt_instance(use_visualization=True) # Use visualized DT to ensure viz manager is ready
⋮----
# Visualize generated surfaces
viz_manager = self._get_viz_manager_instance(output_dir)
surface_files = {}
surface_energies = {}
⋮----
miller_str = ''.join(map(str, analysis_info.miller_index))
slab_file = output_dir / f"slab_{miller_str}.vasp"
# Ensure slab file is written for visualization
⋮----
# Create a dedicated visualization subdirectory
viz_subdir = output_dir / "visualizations"
⋮----
temp_viz_manager = self._get_viz_manager_instance(viz_subdir) # Use temporary manager for subdir
⋮----
viz_paths = temp_viz_manager.visualize_generated_surfaces(
⋮----
# Save SurFF result object
result_path = output_dir / "surff_prediction_result.pkl"
⋮----
"""
        Predicts optimal adsorption sites for a given adsorbate on a surface.

        Args:
            surface_path (str): Absolute path to the surface structure file (e.g., .vasp).
            adsorbate_smi (str): SMILES string or common name for the adsorbate (e.g., "*CO").
            num_sites (int): Number of adsorption sites to sample.
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "02_adsorption").
            llm_review_context (str, optional): Context for LLM if reviewing structures (e.g., "CO oxidation on Pt(111)").

        Returns:
            str: Path to a pickle file containing the AdsorbDiffPredictionOutput object.
                 Also saves the best adsorption configuration to a VASP file in the output directory.
        """
⋮----
dt = self._get_dt_instance(use_visualization=True)
⋮----
# Manually manage initial adsorbate structure, if any, for proper surface expansion
surface_atoms = read(surface_path)
⋮----
# Check if the surface already contains the adsorbate (e.g., from a previous step)
# This is a heuristic to prevent double-adding adsorbate when it's part of the 'surface' input
# For agent 2, input surface should typically be pristine or with other adsorbates, not the one being added.
# However, the underlying dt.predict_adsorption_sites implicitly handles it.
⋮----
adsorb_result = dt.predict_adsorption_sites(
⋮----
# Save best result to VASP
best_config_path = output_dir / "best_adsorption_config.vasp"
⋮----
# Visualize adsorption sites
⋮----
viz_manager = self._get_viz_manager_instance(viz_subdir)
ads_structures = {
# AdsorbDiffOutput contains best_result, so we directly visualize it
⋮----
viz_paths = viz_manager.visualize_adsorption_sites(
⋮----
{ "best": adsorb_result.best_result.energy }, # Use actual energy if available
Path(surface_path).stem # Use filename as identifier
⋮----
# LLM review
⋮----
# LLM review needs the structure and an image
image_for_llm_review = viz_subdir / "llm_review_temp_img.png"
⋮----
llm_review_result = dt._review_structure_with_llm(
⋮----
# Save LLM review result
⋮----
# Save AdsorbDiff result object
result_path = output_dir / "adsorbdiff_prediction_output.pkl"
⋮----
adsorbates_elements_for_mc: List[str], # E.g., ['Pt', 'O'] from the bulk for VSSR-MC to sample
temperature_k: float, # Temperature in Kelvin for MC
⋮----
clean_slab_path: Optional[str] = None, # Optional: path to the clean slab if known, for virtual site generation
surface_indices: Optional[List[int]] = None,  # Slab surface atom indices (free for MC)
adsorbate_indices: Optional[List[int]] = None,  # Pre-existing adsorbate atom indices (free but not "surface")
⋮----
"""
        Simulates surface reconstruction using VSSR-MC for a given surface and adsorbate elements.

        Args:
            surface_with_adsorbate_path (str): Absolute path to the surface structure file (e.g., .vasp)
                                            which may already contain adsorbates.
            adsorbates_elements_for_mc (List[str]): List of element symbols (e.g., ['Pt', 'O'])
                                                    that VSSR-MC can add/remove to simulate reconstruction.
                                                    These are typically the elements in the bulk/surface.
            temperature_k (float): Temperature in Kelvin for the MC simulation.
            total_sweeps (int): Number of Monte Carlo sweeps to perform.
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "03_reconstruction").
            clean_slab_path (str, optional): Absolute path to a clean slab (without adsorbates).
                                             If provided, virtual sites for VSSR-MC will be generated on this.
                                             This helps prevent placing virtual sites above existing adsorbates.

        Returns:
            str: Path to a pickle file containing the VSSRMCResult object.
                 Also saves the lowest energy reconstructed structure to a VASP file.
        """
⋮----
surface_atoms = read(surface_with_adsorbate_path)
⋮----
mc_result = dt.simulate_surface_reconstruction(
⋮----
# Save lowest energy reconstructed structure to VASP
lowest_energy_structure_path = output_dir / "lowest_energy_reconstructed_structure.vasp"
⋮----
# Visualize MC trajectory (if VisualizedGasSolidDigitalTwin is used)
⋮----
# The _run_vssr_mc_subprocess creates trajectory files, we need to collect them
# For simplicity, we just use the raw output of mc_result.structures
# or point to the trajectory if it's saved.
# Here, we generate one from the saved structures.
traj_atoms_list = [s.atoms for s in mc_result.structures]
⋮----
temp_viz_manager = self._get_viz_manager_instance(viz_subdir)
⋮----
# Ensure numpy is imported for linspace
⋮----
# Limit frames for reasonable GIF size and generation time
max_frames = 50
⋮----
indices = np.linspace(0, len(traj_atoms_list) - 1, max_frames, dtype=int)
# Ensure first and last frames are included
if 0 not in indices: indices = np.insert(indices, 0, 0)
if len(traj_atoms_list) - 1 not in indices: indices = np.append(indices, len(traj_atoms_list) - 1)
traj_atoms_list = [traj_atoms_list[i] for i in sorted(list(set(indices)))]
⋮----
mc_gif_path = viz_subdir / f"{Path(surface_with_adsorbate_path).stem}_mc_reconstruction.gif"
⋮----
# Save VSSRMCResult object
result_path = output_dir / "vssr_mc_result.pkl"
⋮----
surface_path: str, # Could be pristine, or with a starting adsorbate
reaction_intermediates: List[str], # E.g., ["*CO", "*O", "*CO2", "*"]
initial_adsorbate_present_on_surface_smi: str, # E.g., "*CO" if surface_path has *CO
⋮----
"""
        Analyzes the reaction pathway by calculating adsorption energies and NEB barriers.

        Args:
            surface_path (str): Absolute path to the reference surface structure file (e.g., lowest energy reconstructed surface).
            reaction_intermediates (List[str]): Ordered list of reaction intermediates (e.g., ["*CO", "*O", "*CO2", "*"]).
            initial_adsorbate_present_on_surface_smi (str): The adsorbate that is already present on the 'surface_path' structure.
                                                          This helps correctly identify surface atoms vs. adsorbate atoms.
                                                          Use "" or None if surface_path is a pristine surface.
            calculate_barriers (bool): Whether to perform NEB calculations for barriers.
            num_sites (int): Number of potential adsorption sites to consider for each intermediate.
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "04_pathway_analysis").

        Returns:
            str: Path to a pickle file containing the CompletePathwayResult object.
                 Also saves the reaction pathway visualization (GIF) and energy diagram (PNG).
        """
⋮----
# Load initial surface atoms to pass to pathway_predictor to help identify surface atoms
initial_surface_atoms = read(surface_path)
⋮----
# Determine current adsorbate indices correctly
# This logic needs to correctly identify atoms belonging to the initial adsorbate from the surface_path
current_adsorbate_indices = []
⋮----
# 🔧 实现吸附物识别：找到最上层的非Pt原子
⋮----
symbols = initial_surface_atoms.get_chemical_symbols()
positions = initial_surface_atoms.positions
⋮----
# 找到所有非Pt原子
non_pt_indices = [i for i, sym in enumerate(symbols) if sym != "Pt"]
⋮----
# 按z坐标排序，取最上层的原子作为吸附物
z_coords = positions[non_pt_indices, 2]
z_threshold = z_coords.max() - 2.0  # 最上层2Å内的原子
⋮----
adsorbate_mask = z_coords > z_threshold
current_adsorbate_indices = [non_pt_indices[i] for i, is_ads in enumerate(adsorbate_mask) if is_ads]
⋮----
pathway_result = dt.analyze_reaction_pathway(
⋮----
surface=initial_surface_atoms, # Pass the Atoms object directly
⋮----
# Visualize reaction pathway (GIF) and energy diagram (PNG)
⋮----
# --- Pathway GIF ---
pathway_structures_for_viz = {}
pathway_adsorption_dir = output_dir / "adsorption" # Pathway predictor saves individual steps here
⋮----
# Collect structures in the correct order for GIF
⋮----
inter_file = pathway_adsorption_dir / f"{inter}_relaxed.vasp"
⋮----
pathway_structures_for_viz[inter] = read(str(inter_file)) # Load Atoms object
⋮----
pathway_gif = viz_manager.visualize_reaction_pathway(
⋮----
# --- Energy Diagram PNG ---
⋮----
# Extract energy profile as done in visualized_digital_twin.py
energies = []
labels = []
is_ts_list = []
⋮----
# The pathway_result's energy_profile attribute directly provides what's needed for visualization
# It's a list of tuples: (label, energy, is_transition_state)
⋮----
energy_diagram_path = viz_manager.visualize_energy_diagram(
⋮----
# Save CompletePathwayResult object
result_path = output_dir / "complete_pathway_result.pkl"
⋮----
pressures: Dict[str, float], # E.g., {"CO": 1.0, "O2": 0.2}
⋮----
"""
        Runs Kinetic Monte Carlo (KMC) simulation based on a previously analyzed reaction pathway.

        Args:
            pathway_result_path (str): Absolute path to a pickle file containing the CompletePathwayResult object.
            temperature_k (float): Temperature in Kelvin for the KMC simulation.
            pressures (Dict[str, float]): Dictionary of gas phase pressures in bar (e.g., {"CO": 1.0, "O2": 0.2}).
            run_id (str): A unique identifier for the current workflow run.
            workflow_step (str): Identifier for the current step (e.g., "05_kmc_simulation").

        Returns:
            str: Path to a pickle file containing the CatMAPResult object.
                 Also saves the KMC dynamics visualization (GIF) if available.
        """
⋮----
pathway_result: CompletePathwayResult = self._load_result_from_pickle(pathway_result_path)
⋮----
kmc_result = dt.run_kmc_simulation(
⋮----
# Visualize KMC dynamics (GIF)
⋮----
# We need the structures for each intermediate to create the GIF
# Fetching the base surface from the pathway_result's initial_atoms if available,
# otherwise try from the reconstructed path.
base_surface_atoms = None
first_step_initial_atoms = None
⋮----
first_step_initial_atoms = getattr(pathway_result.steps[0], "initial_atoms", None)
candidate_base_surfaces = [
⋮----
base_surface_atoms = first_step_initial_atoms
⋮----
base_surface_atoms = read(str(candidate))
⋮----
# Fallback: if no base surface, cannot create dynamic visualization correctly
⋮----
# Get all intermediate structures from pathway analysis output
# This needs to include the base_surface_atoms as the "clean" state if '*' is a state
intermediate_structures = {}
adsorption_candidates = [
pathway_adsorption_dir = None
⋮----
pathway_adsorption_dir = candidate
⋮----
# For the clean surface state, use the base_surface_atoms
⋮----
inter_file = pathway_adsorption_dir / f"{inter}_relaxed.vasp" if pathway_adsorption_dir is not None else None
⋮----
state_history = kmc_result.state_history
time_history = kmc_result.time_history
⋮----
max_frames = 100
total_steps = len(state_history)
⋮----
# Select frames, ensuring first and last are always present
indices = np.linspace(0, total_steps - 1, min(total_steps, max_frames), dtype=int)
⋮----
trajectory = []
titles = []
⋮----
state = state_history[idx]
time = time_history[idx]
⋮----
# Need to copy the Atoms object to avoid modifying it during visualization
atoms_to_add = intermediate_structures[state].copy()
# Optional: Add title to atoms object for visualization_manager
⋮----
kmc_gif_path = viz_subdir / "kmc_surface_dynamics.gif"
⋮----
show_progress=False # Disable verbose progress for tool output
⋮----
# Save CatMAPResult object
result_path = output_dir / "catmap_result.pkl"
⋮----
predictor = self.get_shared_fairchem_predictor(
result = predictor.predict_pathway_energies(
````

## File: camel_agents/tooling/species.py
````python
"""Species parsing, index management, structure serialization, and reaction context normalization."""
⋮----
class WorkflowContextMixin
⋮----
def _set_surface_indices(self, atoms: Atoms, surface_indices: List[int], adsorbate_indices: List[int]) -> None
⋮----
def _update_adsorbate_indices_from_count(self, atoms: Atoms, adsorbate_count: int) -> List[int]
⋮----
adsorbate_indices: List[int] = []
⋮----
adsorbate_indices = list(range(len(atoms) - adsorbate_count, len(atoms)))
surface_indices = [i for i in range(len(atoms)) if i not in adsorbate_indices]
⋮----
"""Infer adsorbate indices after MC reconstruction using reference mapping first.

        Why:
        - Reconstruction can append/reorder atoms (e.g., added slab atoms),
          so tail-based indexing is unsafe.
        - We first map previous adsorbate atoms to reconstructed atoms by
          element + nearest-distance, then fallback to generic heuristics.
        """
expected_count = len(reference_adsorbate_indices or [])
mapped_adsorbate: List[int] = []
⋮----
all_reconstructed_indices = list(range(len(reconstructed_atoms)))
pairs = self._match_adsorbate_atoms_by_element(
⋮----
pair_by_reference: Dict[int, int] = {int(r_idx): int(p_idx) for r_idx, p_idx in pairs}
ordered = [pair_by_reference[r_idx] for r_idx in reference_adsorbate_indices if r_idx in pair_by_reference]
⋮----
seen: set[int] = set()
mapped_adsorbate = []
⋮----
_ = fallback_surface
⋮----
mapped_adsorbate = mapped_adsorbate[:expected_count]
⋮----
surface_indices = [i for i in range(len(reconstructed_atoms)) if i not in set(mapped_adsorbate)]
⋮----
tags = atoms.get_tags() if len(atoms) else np.array([])
⋮----
max_tag = int(np.max(tags))
candidate_ads = [i for i, tag in enumerate(tags.tolist()) if int(tag) == max_tag]
⋮----
surface = [i for i in range(len(atoms)) if i not in set(candidate_ads)]
⋮----
symbols = atoms.get_chemical_symbols()
⋮----
dominant_symbol = Counter(symbols).most_common(1)[0][0]
positions = atoms.get_positions()
⋮----
dominant_indices = [i for i, s in enumerate(symbols) if s == dominant_symbol]
⋮----
surface_z = float(np.max(positions[dominant_indices, 2]))
⋮----
surface_z = float(np.max(positions[:, 2]))
⋮----
adsorbate = [
⋮----
adsorbate = [i for i, s in enumerate(symbols) if s != dominant_symbol]
⋮----
adsorbate = list(range(fallback_surface_count, len(atoms)))
⋮----
surface = [i for i in range(len(atoms)) if i not in set(adsorbate)]
⋮----
clean_atoms = read(state.clean_slab_path)
⋮----
derived_clean = ads_atoms[surface_indices]
derived_path = adsorption_dir / "derived_clean_slab.vasp"
⋮----
def _format_adsorbate_info(self, atoms: Atoms, adsorbate_indices: List[int]) -> str
⋮----
lines = []
⋮----
pos = positions[idx]
⋮----
@staticmethod
    def _read_text_file(path_text: str) -> str
⋮----
path = Path(str(path_text or "").strip())
⋮----
ads_set = set(adsorbate_indices or [])
⋮----
indices = list(range(len(atoms)))
⋮----
indices = [idx for idx in selected_indices if 0 <= idx < len(atoms)]
lines = ["index symbol x y z role"]
⋮----
sym = symbols[idx]
⋮----
role = "adsorbate" if idx in ads_set else "surface"
⋮----
n = max(0, min(int(fixed_count), len(atoms)))
⋮----
digest_lines = [
payload = "\n".join(digest_lines).encode("utf-8")
⋮----
selected = list(range(n))
signature = self._fixed_coordinate_signature(atoms, n)
coord_text = self._format_atom_coordinates_table(
⋮----
def _atoms_to_poscar_text(self, atoms: Atoms) -> str
⋮----
buffer = StringIO()
⋮----
ads_indices = [idx for idx in (adsorbate_indices or []) if 0 <= idx < len(atoms)]
selected = [idx for idx in (selected_indices or []) if 0 <= idx < len(atoms)]
⋮----
serialized_atoms = atoms
header_lines = [f"[{title}]"]
⋮----
serialized_atoms = atoms[selected].copy()
⋮----
selected_syms = [atoms[idx].symbol for idx in selected]
⋮----
poscar_text = self._atoms_to_poscar_text(serialized_atoms)
⋮----
def _build_agent5_structures_payload(self, state: WorkflowState) -> str
⋮----
def _build_clean_surface_from_known_indices(self, base_structure: Atoms, state: WorkflowState) -> Atoms
⋮----
"""Build clean slab using upstream-known index partition (no re-detection)."""
known_surface = list(state.surface_indices or [])
known_ads = list(state.adsorbate_indices or [])
⋮----
clean = base_structure[known_surface].copy()
⋮----
clean = base_structure.copy()
⋮----
valid = [idx for idx in (adsorbate_indices or []) if 0 <= idx < len(base_structure)]
⋮----
positions = base_structure.get_positions()[valid]
anchor_local_idx = int(np.argmin(positions[:, 2]))
anchor = positions[anchor_local_idx]
⋮----
def _species_safe_name(self, species: str) -> str
⋮----
base = self._canonical_species_label(species) or str(species or "").strip()
safe = re.sub(r"[^A-Za-z0-9._-]+", "_", base)
⋮----
"""Agent4 phase-1: tool-based fixed-site product placement for intermediates."""
intermediates = self._get_working_intermediate_sequence(state)
⋮----
placed: Dict[str, Dict[str, Any]] = {}
⋮----
base_copy = base_structure.copy()
⋮----
# If first intermediate is a co-adsorption label, add secondary species
first_label = intermediates[0]
⋮----
primary_ads_indices = list(base_copy.info.get("adsorbate_indices", []))
# Fallback: if adsorbate_indices empty, infer from surface_atom_count
⋮----
primary_ads_indices = list(range(int(state.surface_atom_count), len(base_copy)))
⋮----
site = self._find_nearest_neighbor_site(
⋮----
sec_atoms = self._formula_to_atom_list(sec_formula)
base_pos = np.array(site, dtype=float)
⋮----
angle = 2.0 * np.pi * j / len(sec_atoms)
bond_len = 0.97 if elem == "H" else 1.20
offset_vec = np.array([
⋮----
# Update indices after adding co-adsorbates
surface_count = len(list(state.surface_indices or []))
surface_indices = list(range(surface_count))
adsorbate_indices = list(range(surface_count, len(base_copy)))
⋮----
first_key = self._canonical_species_label(first_label)
⋮----
iter_dir = Path(state.output_base_dir) / state.run_id / "04_pathway" / f"iter_{iteration:02d}_tool_preplacement"
⋮----
# Separate primary species from co-adsorption labels for tool placement.
# E.g. "*CO+*O" → primary "*CO" goes to adsorption tool; "*O" placed later.
tool_intermediates: List[str] = []
coadsorption_map: Dict[str, List[Tuple[str, int]]] = {}  # label → [(formula, count)]
⋮----
# Feed only the primary species to the adsorption tool
primary_label = f"*{primary}"
⋮----
clean_surface = self._build_clean_surface_from_known_indices(base_structure, state)
clean_surface_path = iter_dir / "clean_surface_known_indices.vasp"
⋮----
anchor_site = self._estimate_anchor_from_known_adsorbate(base_structure, list(state.adsorbate_indices or []))
tool_result = self.tools.compute_adsorption_energies(
⋮----
best_configs = tool_result.get("best_configurations", {}) if isinstance(tool_result, dict) else {}
surface_atom_count = len(clean_surface)
surface_template = clean_surface.copy()
⋮----
key = self._canonical_species_label(label)
⋮----
# Determine which tool label was used for the primary species
⋮----
tool_label = f"*{primary}"
⋮----
tool_label = label
⋮----
config = (
⋮----
merged = surface_template.copy()
adsorbate_atoms = config[surface_atom_count:]
⋮----
# Place co-adsorbed secondary species at nearest neighbor sites
⋮----
primary_ads_indices = list(range(surface_atom_count, len(merged)))
⋮----
# Single atom (H, O, etc.) — place directly at site
⋮----
# Multi-atom fragment (OH, H2O, etc.) — build with
# reasonable bond geometry. First atom at site,
# subsequent atoms offset with ~1.0 Å bond length.
⋮----
# Alternate offset directions for each extra atom
⋮----
surface_indices = list(range(surface_atom_count))
adsorbate_indices = list(range(surface_atom_count, len(merged)))
⋮----
product_path = iter_dir / f"{self._species_safe_name(str(label))}_product.vasp"
⋮----
required = [s for s in intermediates if not self._contains_gas_phase_hint(s)]
missing = [s for s in required if self._canonical_species_label(s) not in placed]
⋮----
# Align consecutive intermediates: shift each product's adsorbate so that
# its bottom-most adsorbate atom matches the previous step's bottom-most.
# This ensures the surface-bonded anchor stays at consistent coordinates
# across steps, preventing interpolation collisions in NEB.
⋮----
"""Find nearest-neighbor hollow/bridge site relative to primary adsorbate.

        Uses surface metal atoms near the primary adsorbate centroid to compute
        hollow site positions, then picks the closest one that does NOT overlap
        with any existing adsorbate atom (distance >= *min_dist_from_ads*).
        ``offset`` skips the first N valid sites for placing multiple co-adsorbates.
        """
⋮----
positions = structure.get_positions()
symbols = structure.get_chemical_symbols()
⋮----
# Primary adsorbate centroid
valid_ads = [i for i in primary_ads_indices if i < len(structure)]
⋮----
ads_positions = positions[valid_ads]
centroid_xy = np.mean(ads_positions[:, :2], axis=0)
⋮----
# Surface metal atoms (top layer)
ads_set = set(primary_ads_indices)
surface_indices = [i for i in range(len(structure)) if i not in ads_set]
⋮----
metal_counts = Counter(symbols[i] for i in surface_indices)
dominant = metal_counts.most_common(1)[0][0]
metal_indices = [i for i in surface_indices if symbols[i] == dominant]
⋮----
metal_z = [positions[i][2] for i in metal_indices]
⋮----
top_z = max(metal_z)
top_layer = [i for i in metal_indices if positions[i][2] > top_z - 1.0]
⋮----
# Generate candidate hollow sites from all triplets of nearby top-layer atoms
⋮----
top_positions = np.array([positions[i] for i in top_layer])
candidate_sites: List[Tuple[float, np.ndarray]] = []
⋮----
tri = top_positions[list(combo)]
# Only consider compact triangles (max edge < 3.5 Å, typical fcc nearest-neighbor)
edges = [float(np.linalg.norm(tri[i] - tri[j])) for i, j in [(0, 1), (0, 2), (1, 2)]]
⋮----
site_xy = np.mean(tri[:, :2], axis=0)
site_z = float(np.mean(tri[:, 2])) + adsorption_height
site_pos = np.array([site_xy[0], site_xy[1], site_z])
⋮----
# Distance from adsorbate centroid (for sorting: prefer near but not under)
d_centroid = float(np.linalg.norm(site_xy - centroid_xy))
⋮----
# Minimum distance from ALL adsorbate atoms (must not overlap)
d_min_ads = float(min(np.linalg.norm(site_pos - ads_positions[j]) for j in range(len(ads_positions))))
⋮----
# Sort by distance from centroid (nearest neighbor first)
⋮----
# Skip `offset` sites for placing multiple co-adsorbates at different positions
idx = min(offset, len(candidate_sites) - 1)
⋮----
@staticmethod
    def _formula_to_atom_list(formula: str) -> List[str]
⋮----
"""Parse a chemical formula into a list of element symbols.

        ``"OH"`` → ``["O", "H"]``, ``"H2"`` → ``["H", "H"]``, ``"O"`` → ``["O"]``.
        """
atoms: List[str] = []
⋮----
elem = m.group(1)
count = int(m.group(2)) if m.group(2) else 1
⋮----
"""Align adsorbate positions across consecutive intermediates.

        For each pair of consecutive non-gas intermediates, shift the later
        one's adsorbate atoms so that its bottom-most adsorbate atom aligns
        with the previous intermediate's bottom-most adsorbate atom.  All other
        adsorbate atoms keep their relative coordinates (rigid shift).
        Surface atoms are never modified.

        This fixes the independent-relaxation problem where the same anchor atom
        (e.g. C bonded to surface) ends up at different positions in different
        intermediates.
        """
⋮----
prev_key = None
⋮----
key = self._canonical_species_label(species)
entry = placed.get(key)
⋮----
prev_key = key
⋮----
prev_entry = placed.get(prev_key)
⋮----
"""Shift target adsorbate so its anchor atom aligns with ref's anchor atom.

        Anchor selection priority (general, not reaction-specific):
        1. Heaviest common non-H adsorbate element with lowest z (e.g., C in both CO2 and CO)
        2. Lowest-z adsorbate atom (fallback)
        """
⋮----
ref_struct: Atoms = ref_entry["structure"]
target_struct: Atoms = target_entry["structure"]
ref_ads = list(ref_entry.get("adsorbate_indices", []))
target_ads = list(target_entry.get("adsorbate_indices", []))
⋮----
ref_positions = ref_struct.get_positions()
target_positions = target_struct.get_positions()
ref_symbols = ref_struct.get_chemical_symbols()
target_symbols = target_struct.get_chemical_symbols()
⋮----
ref_valid = [i for i in ref_ads if 0 <= i < len(ref_struct)]
target_valid = [i for i in target_ads if 0 <= i < len(target_struct)]
⋮----
# Find common non-H elements between ref and target adsorbates
ref_heavy = {ref_symbols[i] for i in ref_valid if ref_symbols[i] != "H"}
target_heavy = {target_symbols[i] for i in target_valid if target_symbols[i] != "H"}
common_heavy = ref_heavy & target_heavy
⋮----
# Pick heaviest common element (C > N > O > S by atomic number)
⋮----
anchor_elem = max(common_heavy, key=lambda e: atomic_numbers.get(e, 0))
ref_candidates = [i for i in ref_valid if ref_symbols[i] == anchor_elem]
target_candidates = [i for i in target_valid if target_symbols[i] == anchor_elem]
ref_bottom_idx = min(ref_candidates, key=lambda i: ref_positions[i, 2])
target_bottom_idx = min(target_candidates, key=lambda i: target_positions[i, 2])
⋮----
# Fallback: lowest-z atom
ref_bottom_idx = min(ref_valid, key=lambda i: ref_positions[i, 2])
target_bottom_idx = min(target_valid, key=lambda i: target_positions[i, 2])
⋮----
ref_anchor = ref_positions[ref_bottom_idx]
target_anchor = target_positions[target_bottom_idx]
⋮----
# Compute shift: move target anchor to ref anchor position
shift = ref_anchor - target_anchor
⋮----
# Only shift if displacement is significant (> 0.1 Å)
⋮----
# Apply rigid shift to ALL target adsorbate atoms (surface stays fixed)
new_positions = target_positions.copy()
⋮----
# Update signature since positions changed
⋮----
lines: List[str] = []
⋮----
entry = preplaced_products.get(key)
⋮----
struct = entry.get("structure")
ads_indices = list(entry.get("adsorbate_indices", []))
centroid_text = "(n/a)"
⋮----
centroid = np.mean(struct.get_positions()[ads_indices], axis=0)
centroid_text = f"({centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f})"
⋮----
raw = str(species or "").strip()
⋮----
@staticmethod
    def _parse_coadsorption_label(label: str) -> Tuple[str, List[Tuple[str, int]]]
⋮----
"""Parse co-adsorption labels like ``*CO+*O``, ``*CO2+2*H``.

        Returns (primary_formula, [(secondary_formula, count), ...]).
        For plain labels like ``*CO`` returns (``CO``, []).
        """
raw = str(label).strip().replace(" ", "")
# Split on '+' but only between starred species
parts = re.split(r"\+(?=\d*\*)", raw)
⋮----
# No co-adsorption
cleaned = re.sub(r"^\*+", "", raw)
cleaned = re.sub(r"\([^\)]*\)$", "", cleaned)
⋮----
primary = re.sub(r"^\*+", "", parts[0])
primary = re.sub(r"\([^\)]*\)$", "", primary).upper()
⋮----
secondaries: List[Tuple[str, int]] = []
⋮----
m = re.match(r"(\d*)\*?([A-Za-z0-9]+)", part)
⋮----
count = int(m.group(1)) if m.group(1) else 1
formula = m.group(2).upper()
⋮----
@staticmethod
    def _is_coadsorption_label(label: str) -> bool
⋮----
"""Return True if label contains co-adsorption notation (``+*``)."""
⋮----
@staticmethod
    def _canonical_species_label(label: Optional[str]) -> str
⋮----
raw = str(label).strip()
is_gas_phase = bool(
compact = raw.replace(" ", "")
⋮----
# Preserve co-adsorption labels: *CO+*O → "CO+O", *CO2+2*H → "CO2+2H"
⋮----
parts = re.split(r"\+(?=\d*\*)", compact)
normalized_parts = []
⋮----
cleaned = re.sub(r"\*+", "", part)  # remove all * (including mid-string like 2*H)
⋮----
starred = re.findall(r"\*[A-Za-z0-9]+(?:\*[A-Za-z0-9]+)?", compact)
⋮----
compact = starred[0]
⋮----
compact = raw.split(sep, 1)[0].strip().replace(" ", "")
⋮----
cleaned = re.sub(r"^\*+", "", compact)
⋮----
normalized = cleaned.upper()
⋮----
@staticmethod
    def _extract_primary_adsorbate_species(text: Optional[str]) -> str
⋮----
raw = str(text).strip()
⋮----
head = raw.split(sep, 1)[0].strip()
⋮----
tokens = raw.split()
⋮----
@staticmethod
    def _ensure_star_notation(species: str) -> str
⋮----
label = str(species or "").strip()
⋮----
def _extract_intermediate_sequence_from_text(self, description: str) -> List[str]
⋮----
raw = str(description or "")
⋮----
normalized = raw.replace("→", "->")
chain_candidates = re.findall(r"([*A-Za-z0-9()\-]+(?:\s*->\s*[*A-Za-z0-9()\-]+)+)", normalized)
best: List[str] = []
⋮----
parts = [p.strip(" \t\n,.;") for p in candidate.split("->") if p.strip()]
⋮----
best = parts
⋮----
parts = [p.strip(" \t\n,.;") for p in normalized.split("->") if p.strip()]
⋮----
cleaned: List[str] = []
⋮----
species = self._extract_primary_adsorbate_species(token)
⋮----
deduped: List[str] = []
prev_key = ""
⋮----
key = self._canonical_species_label(item)
⋮----
def _fallback_reaction_context(self, reaction_description: str) -> ReactionContext
⋮----
_ = reaction_description
⋮----
def _fallback_pathway_design(self, state: WorkflowState) -> PathwayDesign
⋮----
_ = state
⋮----
def _normalize_reaction_context(self, state: WorkflowState) -> None
⋮----
context = state.reaction_context
initial = self._extract_primary_adsorbate_species(context.initial_adsorbate)
⋮----
normalized_intermediates: List[str] = []
⋮----
species = self._extract_primary_adsorbate_species(item)
⋮----
initial_key = self._canonical_species_label(initial)
⋮----
normalized_intermediates = [initial]
⋮----
starred_intermediates = [self._ensure_star_notation(item) for item in deduped]
parsed_from_text = self._extract_intermediate_sequence_from_text(state.reaction_description)
⋮----
# User-provided explicit pathway in natural language has highest priority.
starred_intermediates = parsed_from_text
⋮----
def _derive_intermediate_sequence_from_design(self, state: WorkflowState) -> List[str]
⋮----
design = state.pathway_design
⋮----
sequence: List[str] = []
⋮----
reactant = str(spec.reactant_formula or "").strip()
product = str(spec.product_formula or "").strip()
⋮----
def _canonical_step_pair(self, reactant: str, product: str) -> Tuple[str, str]
⋮----
reactant_primary = self._extract_primary_adsorbate_species(reactant)
product_primary = self._extract_primary_adsorbate_species(product)
⋮----
def _formula_element_key(self, label: str) -> Tuple[str, ...]
⋮----
"""Convert a species label to a sorted element tuple for fuzzy matching.

        Both 'CHOH' and 'CH2O' produce ('C', 'H', 'H', 'O').
        """
tokens = self._extract_species_tokens(label)
⋮----
def _align_pathway_design_to_context(self, state: WorkflowState, design: PathwayDesign) -> PathwayDesign
⋮----
context_seq = list(state.reaction_context.intermediates) if state.reaction_context else list(state.intermediates or [])
⋮----
required_pairs = [(context_seq[i], context_seq[i + 1]) for i in range(len(context_seq) - 1)]
⋮----
pair_to_step: Dict[Tuple[str, str], PathwayStepSpec] = {}
# Also build element-count-based fallback index
elem_pair_to_step: Dict[Tuple[Tuple[str, ...], Tuple[str, ...]], PathwayStepSpec] = {}
⋮----
pair = self._canonical_step_pair(spec.reactant_formula, spec.product_formula)
⋮----
elem_pair = (
⋮----
aligned_steps: List[PathwayStepSpec] = []
design_steps = list(design.steps or [])
⋮----
# Try pair-based matching first
all_matched = True
⋮----
expected_pair = self._canonical_step_pair(reactant_label, product_label)
candidate = pair_to_step.get(expected_pair)
⋮----
# Fallback: match by element composition (CHOH == CH2O)
⋮----
expected_elem_pair = (
candidate = elem_pair_to_step.get(expected_elem_pair)
⋮----
all_matched = False
⋮----
# Fallback: sequential matching when LLM uses NEB-endpoint formulas (X->X)
⋮----
_log = logging.getLogger("catdt.species")
⋮----
aligned_steps = []
⋮----
candidate = design_steps[idx]
⋮----
aligned = PathwayDesign(
⋮----
def _infer_formula_element_deltas(self, reactant_formula: str, product_formula: str) -> Tuple[List[str], List[str]]
⋮----
reactant_tokens = Counter(self._extract_species_tokens(reactant_formula))
product_tokens = Counter(self._extract_species_tokens(product_formula))
⋮----
to_add: List[str] = []
to_remove: List[str] = []
⋮----
delta = product_tokens.get(element, 0) - reactant_tokens.get(element, 0)
⋮----
def _count_explicit_remove_elements(self, atoms_to_remove: List[Any]) -> Counter
⋮----
counts: Counter = Counter()
⋮----
explicit_symbol = item.get("element") or item.get("symbol") or item.get("species")
tokens = self._extract_species_tokens(str(explicit_symbol or ""))
⋮----
count_raw = item.get("count", 1)
⋮----
n = max(int(count_raw), 1)
⋮----
n = 1
⋮----
def _normalize_pathway_step_operations(self, design: PathwayDesign) -> PathwayDesign
⋮----
normalized_steps: List[PathwayStepSpec] = []
⋮----
step_copy = PathwayStepSpec(
⋮----
@staticmethod
    def _is_adsorbate_species(label: Optional[str]) -> bool
⋮----
def _get_working_intermediate_sequence(self, state: WorkflowState) -> List[str]
⋮----
context_seq = list(state.reaction_context.intermediates) if state.reaction_context else []
design_seq = self._derive_intermediate_sequence_from_design(state)
⋮----
chosen = context_seq
⋮----
chosen = design_seq
⋮----
chosen = context_seq or design_seq
⋮----
def _resolve_removal_indices(self, step: PathwayStepSpec, adsorbate_indices: List[int]) -> List[int]
⋮----
resolved: List[int] = []
⋮----
idx_candidates: List[Any] = []
⋮----
idx_candidates = [item.get("index")]
⋮----
idx_candidates = list(item.get("indices", []))
⋮----
idx_candidates = [item]
⋮----
resolved_idx = self._resolve_atom_index(
⋮----
explicit_indices = self._resolve_removal_indices(step, adsorbate_indices)
chosen: List[int] = list(explicit_indices)
used = set(chosen)
⋮----
remove_targets = self._count_explicit_remove_elements(step.atoms_to_remove)
⋮----
anchor = np.mean(positions[adsorbate_indices], axis=0) if adsorbate_indices else np.mean(positions, axis=0)
⋮----
existing = sum(1 for idx in chosen if 0 <= idx < len(structure) and structure[idx].symbol == element)
need = max(0, int(target_count) - existing)
⋮----
candidates = [
````

## File: camel_agents/__init__.py
````python
"""CAMEL-based multi-agent package for CatDT.

Keep package import lightweight: delay heavy CAMEL/workflow imports until actually used.
"""
⋮----
__all__ = [
⋮----
# Mechanism search schemas
⋮----
def __getattr__(name)
⋮----
mapping = {
````

## File: camel_agents/adaptive_parameters.py
````python
"""
自适应参数调整器 (AdaptiveParameterTuner)

根据中间结果自动调整计算参数，包括：
1. MC 温度调整
2. NEB 参数优化
3. 收敛判断和参数重试
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
class ConvergenceStatus(Enum)
⋮----
"""收敛状态"""
CONVERGED = "converged"
IMPROVING = "improving"
OSCILLATING = "oscillating"
STUCK = "stuck"
DIVERGING = "diverging"
UNKNOWN = "unknown"
⋮----
@dataclass
class ParameterAdjustment
⋮----
"""参数调整建议"""
parameter: str
current_value: Any
suggested_value: Any
reason: str
confidence: float  # 0-1
⋮----
class MCTemperatureTuner
⋮----
"""
    MC 温度自适应调整器
    
    根据能量轨迹自动调整温度，确保：
    - 足够的探索（能量变化明显）
    - 不过度震荡
    - 合理的接受率（20-60%）
    """
⋮----
# 目标接受率范围
TARGET_ACCEPTANCE_MIN = 0.20
TARGET_ACCEPTANCE_MAX = 0.60
⋮----
# 温度调整因子
TEMP_ADJUST_UP = 1.5
TEMP_ADJUST_DOWN = 0.7
⋮----
adjustment_interval: int = 100  # 每多少步调整一次
⋮----
def update(self, energy: float, accepted: bool, step: int)
⋮----
"""更新状态"""
⋮----
# 按间隔调整温度
⋮----
def _adjust_temperature(self, current_step: int)
⋮----
"""调整温度"""
window_size = min(self.adjustment_interval, len(self.acceptance_history))
⋮----
# 计算最近窗口的接受率
recent_acceptance = sum(self.acceptance_history[-window_size:]) / window_size
⋮----
# 计算能量变化
recent_energies = self.energy_history[-window_size:]
energy_std = np.std(recent_energies) if len(recent_energies) > 1 else 0
⋮----
old_temp = self.temperature
adjustment_reason = ""
⋮----
# 策略 1: 接受率过低 -> 升温
⋮----
adjustment_reason = f"acceptance rate {recent_acceptance:.2%} too low"
⋮----
# 策略 2: 接受率过高 -> 降温
⋮----
adjustment_reason = f"acceptance rate {recent_acceptance:.2%} too high"
⋮----
# 策略 3: 能量无变化（卡住）-> 大幅升温
elif energy_std < 0.001:  # 能量几乎不变
⋮----
adjustment_reason = f"energy stuck (std={energy_std:.6f})"
⋮----
# 策略 4: 能量震荡过大 -> 小幅降温
elif energy_std > 1.0:  # 能量变化太大
⋮----
adjustment_reason = f"energy oscillating (std={energy_std:.3f})"
⋮----
def get_convergence_status(self, window_size: int = 50) -> ConvergenceStatus
⋮----
"""判断收敛状态"""
⋮----
# 检查是否收敛（能量变化很小）
energy_range = max(recent_energies) - min(recent_energies)
if energy_range < 0.01:  # 能量范围 < 10 meV
⋮----
# 检查趋势
early_avg = np.mean(recent_energies[:window_size//2])
late_avg = np.mean(recent_energies[window_size//2:])
⋮----
if late_avg < early_avg - 0.1:  # 能量明显下降
⋮----
if abs(late_avg - early_avg) < 0.01:  # 能量几乎不变
⋮----
# 检查震荡
energy_diffs = np.diff(recent_energies)
sign_changes = sum(1 for i in range(1, len(energy_diffs))
if sign_changes > window_size * 0.6:  # 频繁变号
⋮----
def suggest_restart_parameters(self) -> Dict[str, Any]
⋮----
"""建议重启参数（如果当前参数无效）"""
status = self.get_convergence_status()
⋮----
"sweep_multiplier": 2.0  # 增加步数
⋮----
class NEBParameterOptimizer
⋮----
"""
    NEB 参数优化器
    
    根据 NEB 收敛情况自动调整参数：
    - 图片数量
    - spring constant
    - climbing image 策略
    """
⋮----
def __init__(self)
⋮----
"""
        分析 NEB 失败原因并给出调整建议
        
        Returns:
            参数调整建议列表
        """
adjustments = []
⋮----
# 获取 NEB 结果信息
fmax = getattr(neb_result, 'fmax', None)
n_steps = getattr(neb_result, 'n_steps', 0)
energies = getattr(neb_result, 'energies', [])
⋮----
# 失败类型 1: fmax 爆炸（力太大）
⋮----
current_value=7,  # 默认值
⋮----
suggested_value=0.05,  # 降低弹簧常数
⋮----
# 失败类型 2: 不收敛（步数用尽）
elif n_steps >= 300:  # 假设最大 300 步
⋮----
suggested_value=0.1,  # 放宽收敛标准
⋮----
# 失败类型 3: 能量不单调（可能错过了 TS）
⋮----
# 检查是否有多个峰值
peaks = sum(1 for i in range(1, len(energies)-1)
⋮----
# 检查是否存在能量阱（不合理）
valleys = sum(1 for i in range(1, len(energies)-1)
⋮----
suggested_value="idpp",  # Image Dependent Pair Potential
⋮----
# 失败类型 4: 初末态结构问题（从验证器获取）
# 这部分在 enhanced_neb_validator 中处理
⋮----
"""
        根据失败历史获取优化后的参数
        
        Parameters:
            attempt_number: 当前尝试次数
            previous_failures: 之前失败的记录
        
        Returns:
            优化后的参数字典
        """
params = {
⋮----
# 根据尝试次数逐步调整
⋮----
# 第一次：标准参数
⋮----
# 第二次：增加图片数，使用更好的插值
⋮----
# 第三次：更多图片，CI-NEB，更宽松的收敛
⋮----
# 更多尝试：激进调整
⋮----
# 根据历史失败模式调整
⋮----
adjustments = self.analyze_failure(
⋮----
class AdaptiveNEBRunner
⋮----
"""
    自适应 NEB 运行器
    
    自动尝试不同的参数组合直到成功
    """
⋮----
def __init__(self, max_attempts: int = 5)
⋮----
"""
        运行 NEB，自动重试不同的参数
        
        Parameters:
            neb_fn: NEB 计算函数 (params) -> result
            initial: 初态结构
            final: 末态结构
            validate_fn: 结果验证函数
        
        Returns:
            (result, metadata)
        """
⋮----
# 获取优化参数
params = self.optimizer.get_optimized_parameters(
⋮----
# 运行 NEB
result = neb_fn(params)
⋮----
# 验证结果
⋮----
# 成功！
metadata = {
⋮----
# 所有尝试都失败
⋮----
# =============================================================================
# 集成到工作流
⋮----
class AdaptiveWorkflowParameters
⋮----
"""
    自适应工作流参数管理
    
    根据前面的步骤结果自动调整后续参数
    """
⋮----
def __init__(self, base_config: Dict[str, Any])
⋮----
"""
        根据表面复杂度调整参数
        
        复杂表面（大晶胞、多元素）需要更多 MC 步数
        """
config = self.base_config.copy()
⋮----
# 获取表面信息
n_atoms = getattr(surface_result, 'n_atoms', 0)
n_elements = len(set(getattr(surface_result, 'symbols', [])))
⋮----
# 调整 MC 步数
⋮----
"""
        根据反应复杂度调整参数
        
        复杂反应（多步骤、高能量变化）需要更精细的 NEB
        """
⋮----
n_steps = len(getattr(pathway_result, 'steps', []))
⋮----
# 多步骤反应，减少每步的 NEB 图片数以节省时间
config['neb_n_images'] = 5  # 而不是默认 7
⋮----
def get_adjustment_summary(self) -> str
⋮----
"""获取调整摘要"""
⋮----
lines = ["Parameter Adjustments:"]
⋮----
# 便捷函数
⋮----
"""
    便捷的 MC 参数自动调优
    
    Example:
        >>> energies = [...]  # MC 能量历史
        >>> accepted = [...]  # 每一步是否接受
        >>> new_params = auto_tune_mc_parameters(energies, accepted, 500.0)
        >>> print(new_params['temperature'])
    """
tuner = MCTemperatureTuner(initial_temperature=current_temp)
⋮----
status = tuner.get_convergence_status()
restart_params = tuner.suggest_restart_parameters()
````

## File: camel_agents/camel_model_backend.py
````python
#!/usr/bin/env python
"""CAMEL model backend utilities."""
⋮----
def _parse_custom_token_limits() -> Dict[str, int]
⋮----
raw = str(os.getenv("CATDT_CUSTOM_MODEL_TOKEN_LIMITS", "")).strip()
⋮----
obj = json.loads(raw)
⋮----
out: Dict[str, int] = {}
⋮----
def _patch_unified_model_token_limit() -> None
⋮----
"""
    Patch CAMEL UnifiedModelType.token_limit for custom model aliases.

    Why:
    - CAMEL treats unknown model strings as UnifiedModelType and logs:
      "Unknown model ... context window size not defined".
    - We use OPENAI-compatible model names (e.g. `gpt-5.4`, `claude-opus-4-6`).
    - Add a deterministic local mapping to avoid noisy warning and undefined limit.
    """
⋮----
original_prop = getattr(UnifiedModelType, "token_limit", None)
original_getter = getattr(original_prop, "fget", None)
⋮----
default_map: Dict[str, int] = {
⋮----
def _patched_token_limit(self: Any) -> int
⋮----
model_name = str(self)
⋮----
UnifiedModelType._catdt_token_limit_patched = True  # type: ignore[attr-defined]
UnifiedModelType._catdt_original_token_limit_getter = original_getter  # type: ignore[attr-defined]
UnifiedModelType.token_limit = property(_patched_token_limit)  # type: ignore[assignment]
⋮----
def _build_model_config(temperature: float, **kwargs: Any) -> Dict[str, Any]
⋮----
model_config: Dict[str, Any] = {"temperature": temperature}
⋮----
def _build_token_counter(model_name: str)
⋮----
"""Build a token counter that tracks the configured Claude model itself.

    CAMEL's `OpenAICompatibleModel.token_counter` currently hardcodes
    `OpenAITokenCounter(ModelType.GPT_4O_MINI)` when no counter is injected.
    For our strict `claude-opus-4-6` backend this is architecturally wrong and
    causes the GPT tokenizer path to appear in runtime diagnostics. Inject a
    model-aware counter explicitly at backend creation time.
    """
⋮----
"""Create a strict CAMEL model backend (no fallback)."""
⋮----
model_name = str(model or os.getenv("OPENAI_MODEL", "gpt-5.4") or "").strip()
⋮----
key = api_key or os.getenv("OPENAI_API_KEY")
url = base_url or os.getenv("OPENAI_BASE_URL")
⋮----
platform = ModelPlatformType.OPENAI_COMPATIBLE_MODEL if url else ModelPlatformType.DEFAULT
request_timeout = kwargs.pop("timeout", None)
⋮----
request_timeout = float(os.getenv("OPENAI_REQUEST_TIMEOUT_SEC", "1800"))
⋮----
max_retries = kwargs.pop("max_retries", None)
⋮----
max_retries = int(os.getenv("OPENAI_MAX_RETRIES", "3"))
⋮----
env_seed = str(os.getenv("OPENAI_SEED", "")).strip()
````

## File: camel_agents/cross_run_cache.py
````python
"""
跨运行结果缓存系统 (CrossRunCache)

避免重复计算相同结构的 SurFF、AdsorbDiff 等耗时步骤。
基于结构指纹和参数哈希实现智能缓存。
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
class StructureFingerprinter
⋮----
"""结构指纹生成器 - 生成结构的唯一标识"""
⋮----
@staticmethod
    def generate_fingerprint(atoms: Atoms, tolerance: float = 0.01) -> str
⋮----
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
⋮----
# 按元素和z坐标排序
sort_indices = np.lexsort((positions[:, 2], positions[:, 1], positions[:, 0], symbols))
sorted_symbols = symbols[sort_indices]
sorted_positions = positions[sort_indices]
⋮----
# 量化位置（容差）
quantized_positions = np.round(sorted_positions / tolerance).astype(int)
⋮----
# 构建指纹数据
fingerprint_data = {
⋮----
# 计算最终指纹
fp_string = json.dumps(fingerprint_data, sort_keys=True)
⋮----
@staticmethod
    def generate_params_fingerprint(params: Dict) -> str
⋮----
"""生成参数字典的指纹"""
# 过滤掉不可序列化的参数
serializable = {}
⋮----
param_string = json.dumps(serializable, sort_keys=True)
⋮----
@dataclass
class CacheEntry
⋮----
"""缓存条目"""
cache_key: str
operation: str  # "surff", "adsorbdiff", "neb", etc.
structure_fingerprint: str
params_fingerprint: str
result_path: str
timestamp: datetime
access_count: int = 0
last_accessed: Optional[datetime] = None
⋮----
def to_dict(self) -> Dict
⋮----
@classmethod
    def from_dict(cls, data: Dict) -> "CacheEntry"
⋮----
class CrossRunCache
⋮----
"""
    跨运行结果缓存
    
    支持 SQLite 后端，持久化存储计算结果。
    自动管理缓存大小和过期时间。
    """
⋮----
# 初始化数据库
⋮----
# 统计
⋮----
def _init_db(self)
⋮----
"""初始化 SQLite 数据库"""
⋮----
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
⋮----
# 查询数据库
⋮----
cursor = conn.execute(
row = cursor.fetchone()
⋮----
timestamp = datetime.fromisoformat(timestamp_str)
⋮----
# 检查过期
ttl = ttl_days or self.default_ttl_days
⋮----
# 加载结果
⋮----
result = self._load_result(result_path)
⋮----
# 更新访问统计
⋮----
"""
        保存结果到缓存
        
        Returns:
            缓存键
        """
⋮----
# 保存结果文件
result_path = self._save_result(cache_key, result)
⋮----
# 保存到数据库
⋮----
# 检查缓存大小
⋮----
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
⋮----
cached = self.get(operation, structure, params, ttl_days)
⋮----
# 执行计算
⋮----
start_time = datetime.now()
result = compute_fn()
compute_time = (datetime.now() - start_time).total_seconds()
⋮----
# 保存到缓存
⋮----
# 记录节省的时间（下次使用）
⋮----
"""
        使缓存失效
        
        - operation=None, structure=None: 清空所有缓存
        - operation="surff", structure=None: 清空所有 surff 缓存
        - operation="surff", structure=bulk: 清空该结构的 surff 缓存
        """
⋮----
# 清空所有
⋮----
# 清空特定操作
⋮----
# 清空特定结构和操作
⋮----
def get_stats(self) -> Dict
⋮----
"""获取缓存统计"""
⋮----
cursor = conn.execute("SELECT COUNT(*) FROM cache_entries")
total_entries = cursor.fetchone()[0]
⋮----
operation_counts = {row[0]: row[1] for row in cursor.fetchall()}
⋮----
hit_rate = 0
total_requests = self.stats["hits"] + self.stats["misses"]
⋮----
hit_rate = self.stats["hits"] / total_requests
⋮----
def _save_result(self, cache_key: str, result: Any) -> str
⋮----
"""保存结果到文件"""
# 使用两级目录结构避免单个目录文件过多
subdir = cache_key[:2]
result_dir = self.cache_dir / "results" / subdir
⋮----
result_path = result_dir / f"{cache_key}.pkl"
⋮----
def _load_result(self, result_path: str) -> Any
⋮----
"""从文件加载结果"""
⋮----
def _remove_entry(self, cache_key: str)
⋮----
"""删除缓存条目"""
⋮----
result_path = Path(row[0])
⋮----
def _enforce_size_limit(self)
⋮----
"""强制执行大小限制"""
current_size = self._get_cache_size_mb()
max_size_mb = self.max_size_gb * 1024
⋮----
# 删除最旧的条目
⋮----
while current_size > max_size_mb * 0.9:  # 清理到 90%
⋮----
def _get_cache_size_mb(self) -> float
⋮----
"""获取缓存大小（MB）"""
total_size = 0
⋮----
# Forward declaration to avoid circular import
class CachedCatDTTools
⋮----
"""
    带缓存的 CatDTTools 包装器
    
    自动缓存 SurFF、AdsorbDiff 等耗时操作的结果。
    """
⋮----
def __init__(self, tools: Any, cache: Optional[CrossRunCache] = None)
⋮----
"""缓存版本的 generate_surfaces"""
⋮----
bulk = read(bulk_structure_path)
params = {"top_n": top_n_surfaces, **kwargs}
⋮----
def compute()
⋮----
result_path = self.cache.compute_or_cache(
⋮----
# 复制结果到当前运行目录
⋮----
current_output = Path(self.tools.output_base_dir) / run_id / workflow_step
⋮----
# 读取缓存的结果并重新保存到当前位置
result = self.tools._load_result_from_pickle(result_path)
new_path = current_output / "surff_prediction_result.pkl"
⋮----
"""缓存版本的 predict_adsorption_sites"""
⋮----
surface = read(surface_path)
params = {"adsorbate": adsorbate_smi, "num_sites": num_sites, **kwargs}
⋮----
# 复制到当前位置
⋮----
run_id = kwargs.get('run_id', 'default_run')
workflow_step = kwargs.get('workflow_step', '02_adsorption')
⋮----
new_path = current_output / "adsorbdiff_prediction_result.pkl"
⋮----
def get_cache_stats(self) -> Dict
⋮----
# =============================================================================
# 便捷函数
⋮----
def clear_all_cache(cache_dir: str = ".catdt_cache")
⋮----
"""清空所有缓存"""
cache = CrossRunCache(cache_dir=cache_dir)
⋮----
def show_cache_stats(cache_dir: str = ".catdt_cache")
⋮----
"""显示缓存统计"""
⋮----
stats = cache.get_stats()
````

## File: camel_agents/gas_solid_digital_twin.py
````python
"""
Gas-Solid Catalysis Digital Twin - 气固界面催化数字孪生

这个模块提供了完整的气固界面催化反应模拟流程，从bulk结构到最终的KMC模拟。

完整流程：
1. 从bulk结构生成最可能暴露的表面 (SurFF)
2. 在每个表面上预测初始反应物的吸附位点 (AdsorbDiff)
3. 模拟表面重构 (VSSR-MC, 不考虑电势)
4. 计算反应路径的自由能和能垒 (Pathway/Fairchem)
5. 进行动力学蒙特卡洛模拟 (KMC/CatMAP)

使用方式:
    from camel_agents.gas_solid_digital_twin import GasSolidDigitalTwin

    # 初始化
    dt = GasSolidDigitalTwin(
        surff_root="deps/SurFF",
        adsorbdiff_root="deps/AdsorbDiff",
        surface_sampling_root="deps/surface-sampling",
        fairchem_root="deps/fairchem",
    )

    # 运行完整流程
    result = dt.run_complete_workflow(
        bulk_structure="POSCAR",
        reaction_intermediates=["*O", "*OH", "*OOH", "*"],
        initial_reactant="*O",
        temperature=500,  # K
        output_dir="output/gas_solid",
    )

    print(result.summary())

作者: Claude
日期: 2026-01-20
"""
⋮----
# 导入各个模块的预测器
⋮----
# =============================================================================
# 数据结构定义
⋮----
@dataclass
class SurfaceAnalysisResult
⋮----
"""单个表面的完整分析结果"""
miller_index: Tuple[int, int, int]
surface_energy: float  # eV/Å²
area_fraction: float  # 表面暴露比例
pristine_surface: Atoms  # 初始表面
adsorption_sites: AdsorbDiffOutput  # 吸附位点预测结果
reconstructed_surface: Optional[VSSRMCResult] = None  # 重构后的表面
pathway_analysis: Optional[CompletePathwayResult] = None  # 反应路径分析
⋮----
def __repr__(self)
⋮----
@dataclass
class GasSolidWorkflowResult
⋮----
"""气固界面催化数字孪生的完整结果"""
bulk_formula: str
reaction_intermediates: List[str]
temperature: float  # K
surff_result: SurFFResult  # SurFF表面预测结果
surface_analyses: List[SurfaceAnalysisResult]  # 每个表面的分析结果
best_surface: Optional[SurfaceAnalysisResult] = None  # 最优表面
kmc_result: Optional[Any] = None  # KMC模拟结果
output_dir: Optional[str] = None
⋮----
def summary(self, top_n: int = 3) -> str
⋮----
"""生成结果摘要"""
lines = [
⋮----
best_site = surface.adsorption_sites.best_result
⋮----
recon = surface.reconstructed_surface.lowest_energy_structure
⋮----
pathway = surface.pathway_analysis
⋮----
rds = pathway.rate_determining_step
⋮----
# 主类
⋮----
class GasSolidDigitalTwin
⋮----
"""
    气固界面催化数字孪生

    整合了从bulk结构到KMC模拟的完整工作流程。
    """
⋮----
"""
        初始化气固界面催化数字孪生

        参数:
            surff_root: SurFF库的路径
            adsorbdiff_root: AdsorbDiff库的路径
            surface_sampling_root: surface-sampling库的路径
            fairchem_root: Fairchem库的路径
            surff_checkpoint: SurFF模型检查点路径（可选）
            adsorbdiff_checkpoint: AdsorbDiff模型检查点路径（可选）
            fairchem_model: Fairchem模型名称 ("uma-s-1p1" 或 "uma-m-1p1")
            fairchem_model_path: Fairchem模型本地路径（可选，如果提供则不从HF下载）
            vssr_mc_model: VSSR-MC使用的力场模型 ("CHGNetNFF" 或 "MACENFF")
            use_gpu: 是否使用GPU
            logger: 日志记录器
            use_llm_controller: 是否使用LLM控制NEB结构准备
            llm_model: LLM模型名称（如"claude-opus-4-5-20251101"）
        """
⋮----
# 初始化各个预测器
⋮----
"""
        步骤1: 从bulk结构生成最可能暴露的表面

        参数:
            bulk_structure: bulk结构文件路径或ASE Atoms对象
            top_n: 生成前N个最可能暴露的表面
            output_dir: 输出目录

        返回:
            (SurFF预测结果, 表面slab列表)
        """
⋮----
# 使用SurFF预测表面
surff_result = self.surff_predictor.predict(
⋮----
# 生成表面slab结构
# 注意: SurFF只预测表面能，需要手动生成slab结构
# 这里我们使用ASE的surface模块生成标准slab
⋮----
bulk = read(bulk_structure)
⋮----
bulk = bulk_structure
⋮----
slabs = []
⋮----
miller = surface.miller_index
slab = self._generate_slab_from_bulk(bulk, miller)
⋮----
# 保存slab结构
⋮----
slab_file = os.path.join(output_dir, f"slab_{miller[0]}{miller[1]}{miller[2]}.vasp")
⋮----
"""
        步骤2: 在表面上预测初始反应物的吸附位点

        如果吸附后分子间距离过近（考虑周期性边界），自动扩大表面并重新吸附

        参数:
            surface: 表面slab结构
            adsorbate: 吸附物（SMILES格式或名称，如"*CO", "*O"）
            num_sites: 尝试的吸附位点数
            output_dir: 输出目录
            min_adsorbate_distance: 吸附原子最小允许距离 (Å)，考虑周期性
            max_expansion_attempts: 最大扩胞尝试次数

        返回:
            吸附位点预测结果
        """
⋮----
current_surface = surface.copy()
n_surface_atoms = len(surface)  # 记录原始表面原子数
expansion_factor = 1  # 追踪扩胞倍数
⋮----
# 进行吸附预测
result = self.adsorbdiff_predictor.predict(
⋮----
# 检查最佳吸附构型中吸附原子的周期性距离
best_structure = result.best_result.final_structure
⋮----
# 扩大表面并重新吸附
⋮----
P = np.diag([expansion_factor, expansion_factor, 1])
current_surface = make_supercell(surface, P)
n_surface_atoms = len(current_surface)  # 更新表面原子数
⋮----
# 距离合适，返回结果
⋮----
"""
        检查吸附分子与其周期性镜像之间的距离

        通用逻辑：
        1. 对比吸附前后识别吸附原子（新增的原子）
        2. 计算吸附分子质心
        3. 检查cell大小，确保周期性镜像不会太近

        Parameters
        ----------
        surface_before : Atoms
            吸附前的表面
        surface_after : Atoms
            吸附后的结构
        threshold : float
            最小允许距离 (Å)

        Returns
        -------
        (is_too_close, min_distance)
        """
n_surface = len(surface_before)
n_total = len(surface_after)
⋮----
# 吸附的原子索引（新增的原子）
ads_indices = list(range(n_surface, n_total))
⋮----
# 获取吸附原子位置
ads_pos = surface_after.get_positions()[ads_indices]
cell = surface_after.cell
⋮----
# 对于吸附分子（无论单原子还是多原子），检查与周期性镜像的距离
# 计算吸附分子的质心位置
ads_center = ads_pos.mean(axis=0)
⋮----
# 最小周期性距离约为 min(cell_x, cell_y) / 2
# 这是分子质心与其周期性镜像之间的最小可能距离
cell_lengths = cell.lengths()
min_cell_length = min(cell_lengths[0], cell_lengths[1])  # xy平面
min_periodic_dist = min_cell_length / 2.0
⋮----
# Debug: 记录cell信息
⋮----
is_too_close = min_periodic_dist < threshold
⋮----
temperature: float = 500.0,  # K
⋮----
"""
        步骤3: 使用VSSR-MC模拟表面重构（不考虑电势）

        参数:
            surface: 表面slab结构
            adsorbates: 可能的吸附物种列表（元素符号，如["O", "H"]）
            temperature: 温度 (K)
            total_sweeps: MC总扫描次数
            output_dir: 输出目录

        返回:
            VSSR-MC模拟结果
        """
⋮----
# Convert temperature from K to kT units for MC sampling
k_B = 8.617333262e-5  # eV/K
mc_temperature = k_B * temperature / 0.025  # Normalize to ~1.0 at 300K
mc_temperature = max(0.5, min(2.0, mc_temperature))
⋮----
# Run VSSR-MC in a subprocess to avoid mpi4py/UCX state corruption
# from other models (SurFF, AdsorbDiff) loaded in this process
result = self._run_vssr_mc_subprocess(
⋮----
"""
        Run VSSR-MC in a subprocess to isolate it from GPU/MPI state corruption.
        """
⋮----
# Create temporary files for communication
⋮----
surface_path = f.name
⋮----
# Use externally provided adsorbate_indices if available,
# otherwise fall back to element-based heuristic detection
⋮----
detected_ads_indices = list(adsorbate_indices)
⋮----
surface_elements = set(self.catalyst_elements) if hasattr(self, 'catalyst_elements') else set()
⋮----
element_counts = Counter(surface.get_chemical_symbols())
surface_elements = {element_counts.most_common(1)[0][0]}
detected_ads_indices = [
⋮----
num_adsorbate_atoms = len(detected_ads_indices)
⋮----
# 创建clean_slab（不含吸附物），用于生成虚拟位点
clean_slab = surface.copy()
⋮----
# 保存clean_slab
⋮----
clean_slab_path = f.name
⋮----
config = {
⋮----
"canonical": True,  # Use canonical mode with fixed adsorbate count
⋮----
# Pass surface/adsorbate indices for proper constraint handling
⋮----
# Electrochemical parameters (liquid-solid Pourbaix mode)
⋮----
config_path = f.name
⋮----
output_path = f.name
⋮----
# Run subprocess
current_dir = Path(__file__).resolve().parent
project_root = current_dir.parent
candidates = [
script_path = next((str(path) for path in candidates if path.exists()), None)
⋮----
candidate_text = ", ".join(str(path) for path in candidates)
⋮----
proc = subprocess.run(
⋮----
# Log subprocess output
⋮----
# Load result
⋮----
output_data = pickle.load(f)
⋮----
error = output_data.get("error", "Unknown error")
tb = output_data.get("traceback", "")
⋮----
# Reconstruct VSSRMCResult
⋮----
def make_structure(s)
⋮----
structures = [make_structure(s) for s in output_data["structures"]]
best_structure = make_structure(output_data["best_structure"])
lowest_structure = make_structure(output_data["lowest_energy_structure"])
⋮----
# Cleanup temporary files
⋮----
"""
        步骤4: 计算反应路径的自由能和能垒

        参数:
            surface: 表面slab结构（可以是重构后的，可能包含吸附物）
            intermediates: 反应中间体列表（按反应顺序，如["*O", "*OH", "*OOH", "*"]）
            calculate_barriers: 是否计算能垒（使用NEB）
            num_sites: 每个中间体尝试的吸附位点数
            n_frames: NEB计算的帧数
            fmax: NEB力收敛标准 (eV/Å)
            output_dir: 输出目录
            current_adsorbate_indices: 当前吸附物的原子索引列表（如果表面已有吸附物）

        返回:
            完整的反应路径分析结果
        """
⋮----
result = self.pathway_predictor.predict_pathway(
⋮----
"""
        步骤5: 进行动力学蒙特卡洛模拟

        参数:
            pathway_result: 反应路径分析结果
            temperature: 温度 (K)
            pressures: 气相物种分压 (bar)，如{"CO": 1.0, "O2": 0.21}
            time_end: 模拟时间 (s)
            output_dir: 输出目录

        返回:
            KMC模拟结果
        """
⋮----
# 导入CatMAP预测器
⋮----
# 创建CatMAP预测器
predictor = CatMAPPredictor(work_dir=output_dir or tempfile.mkdtemp())
⋮----
# 从pathway_result构建反应机理
⋮----
# 当机理中包含元素气相（如 H_g/O_g）时，使用 frozen_gas 避免
# 对不可用分子热化学参数的依赖。
gas_names = list(getattr(predictor, "gases", {}).keys())
⋮----
# 设置条件
⋮----
# 运行单点计算（针对当前表面）
surface_name = pathway_result.surface_formula
result = predictor.run_single_point(surface=surface_name, output_dir=output_dir)
⋮----
"""
        从pathway结果构建CatMAP反应机理

        参数:
            predictor: CatMAPPredictor实例
            pathway_result: 反应路径分析结果
            pressures: 气相物种分压
        """
⋮----
def get_composition(formula: str) -> Counter
⋮----
"""从化学式获取元素组成"""
⋮----
# 处理带*号的吸附物种名称
formula = formula.replace('*', '')
⋮----
# 解析化学式
pattern = r'([A-Z][a-z]?)(\d*)'
composition = Counter()
⋮----
element = match.group(1)
count = int(match.group(2)) if match.group(2) else 1
⋮----
def get_element_change(reactant_formula: str, product_formula: str) -> dict
⋮----
"""计算从反应物到产物的元素变化"""
reactant_comp = get_composition(reactant_formula)
product_comp = get_composition(product_formula)
change = {}
all_elements = set(reactant_comp.keys()) | set(product_comp.keys())
⋮----
diff = product_comp[elem] - reactant_comp[elem]
⋮----
# 定义常见气相物种组成
gas_compositions = {
⋮----
def identify_intact_surface_exchange(reactant_clean: str, product_clean: str) -> Optional[Tuple[str, str]]
⋮----
"""Prefer explicit molecular adsorption/desorption when the adsorbate equals a known gas."""
⋮----
# 收集路径元素：active 表示该元素在步骤间发生了净增减，
# conserved 表示路径中存在但总是守恒。
pathway_elements = set()
active_elements = set()
⋮----
reactant_clean = step.reactant_adsorbate.replace('*', '')
product_clean = step.product_adsorbate.replace('*', '')
⋮----
conserved_elements = set(pathway_elements) - set(active_elements)
⋮----
# 由用户输入分压估计元素活度
element_activity: Dict[str, float] = {}
⋮----
comp = gas_compositions.get(gas_name, get_composition(gas_name))
⋮----
p = max(float(pressure), 1e-20)
⋮----
n = max(int(count), 1)
activity = p ** (1.0 / float(n))
⋮----
# active 元素用于真实计量修正；conserved 元素仅作为“旁观气体”参考，
# 避免 CatMAP 在缺失元素参考时报错。
preferred_active_sources = {
preferred_conserved_sources = {
⋮----
atomic_sources: Dict[str, str] = {}
conserved_reference_gases: Dict[str, str] = {}
effective_pressures: Dict[str, float] = dict(pressures)
⋮----
source = preferred_active_sources.get(elem, elem)
⋮----
activity = max(float(element_activity.get(elem, 1e-12)), 1e-12)
atoms_per_mol = int(gas_compositions[source].get(elem, 1) or 1)
source_pressure = activity ** float(atoms_per_mol)
⋮----
prev_pressure = max(float(effective_pressures.get(source, 0.0)), 0.0)
⋮----
source = preferred_conserved_sources.get(elem, elem)
⋮----
direct_exchange = identify_intact_surface_exchange(reactant_clean, product_clean)
⋮----
gas_comp = get_composition(gas_name)
⋮----
# 添加气相物种（包含用户输入分压和为质量守恒/元素参考补齐的库气体）
⋮----
# 构建反应表达式（包含气相物种）
⋮----
reactant = step.reactant_adsorbate
product = step.product_adsorbate
⋮----
# 分析元素变化
elem_change = get_element_change(reactant, product)
⋮----
# 构建吸附物种表达式
reactant_clean = reactant.replace('*', '')
product_clean = product.replace('*', '')
⋮----
reactant_expr = '*_s'
⋮----
reactant_expr = f"{reactant_clean}*"
⋮----
product_expr = '*_s'
⋮----
product_expr = f"{product_clean}*"
⋮----
gas_species = f"{gas_name}_g"
⋮----
reaction_expr = f"{reactant_expr} + {gas_species} -> {product_expr}"
⋮----
reaction_expr = f"{reactant_expr} -> {product_expr} + {gas_species}"
⋮----
# 构建完整反应表达式，包括气相物种
# 注意：CatMAP 反应表达式中的化学计量系数应为整数。
# 这里使用有理数累积后统一放大，避免出现 0.5H2_g 这类分数字符串。
reactant_stoich: Dict[str, Fraction] = {reactant_expr: Fraction(1, 1)}
product_stoich: Dict[str, Fraction] = {product_expr: Fraction(1, 1)}
⋮----
def add_stoich(target: Dict[str, Fraction], species: str, coeff: Fraction)
⋮----
source_gas = atomic_sources[elem]
⋮----
atoms_per_mol = int(gas_compositions[source_gas].get(elem, 1) or 1)
coeff = Fraction(abs(int(change)), atoms_per_mol)
⋮----
gas_species = f"{source_gas}_g"
if change > 0:  # 需要添加原子（从气相获取）
⋮----
elif change < 0:  # 需要移除原子（释放到气相）
⋮----
denominators = [
scale_factor = 1
⋮----
scale_factor = lcm(scale_factor, int(denom))
⋮----
def render_side(stoich_map: Dict[str, Fraction]) -> List[str]
⋮----
parts: List[str] = []
⋮----
scaled = frac * scale_factor
⋮----
# 防御性分支：避免非整数系数泄漏到 CatMAP 输入。
scaled = Fraction(int(round(float(scaled))), 1)
coeff_int = int(scaled.numerator)
⋮----
reactant_parts = render_side(reactant_stoich)
product_parts = render_side(product_stoich)
⋮----
# 对守恒元素补充“旁观气体”到反应两侧，保证 CatMAP 能构建元素参考。
⋮----
gas_species = f"{source}_g"
⋮----
# 组合反应表达式
reaction_expr = " + ".join(reactant_parts) + " -> " + " + ".join(product_parts)
⋮----
# 添加吸附物种
⋮----
clean_name = ads_name.replace('*', '')
⋮----
# 添加过渡态（如果有能垒信息）
⋮----
ts_name = f"{step.name}-TS"
ts_energy = step.reactant_energy + step.activation_energy
⋮----
# 由 CatMAP 自动从 gas species 推断原子参考。
# 显式 set_atomic_reservoirs 在不同 CatMAP 版本间存在兼容性差异，
# 这里避免手动覆盖以提升稳定性。
⋮----
"""
        运行完整的气固界面催化数字孪生流程

        参数:
            bulk_structure: bulk结构文件路径或ASE Atoms对象
            reaction_intermediates: 反应中间体列表（如["*O", "*OH", "*OOH", "*"]）
            initial_reactant: 初始反应物（如"*O"）
            temperature: 反应温度 (K)
            top_n_surfaces: 分析前N个最可能暴露的表面
            num_adsorption_sites: 每个表面尝试的吸附位点数
            reconstruction_sweeps: VSSR-MC扫描次数
            calculate_barriers: 是否计算反应能垒
            neb_frames: NEB计算的帧数
            neb_fmax: NEB力收敛标准 (eV/Å)
            run_kmc: 是否运行KMC模拟
            pressures: 气相物种分压（如果运行KMC）
            output_dir: 输出目录

        返回:
            完整的工作流程结果
        """
⋮----
# 创建输出目录
⋮----
# 获取bulk结构的化学式
⋮----
bulk_formula = bulk.get_chemical_formula()
⋮----
# 步骤1: 生成表面
⋮----
# For VSSR-MC surface reconstruction, we should only add/remove surface atoms
# NOT the adsorbate molecule atoms (which would break the molecules)
# Extract surface elements from the bulk structure
surface_elements = list(set(bulk.get_chemical_symbols()))
⋮----
# 对每个表面进行完整分析
surface_analyses = []
⋮----
# 创建表面专用输出目录
surface_dir = os.path.join(output_dir, f"surface_{surface_info.miller_index[0]}{surface_info.miller_index[1]}{surface_info.miller_index[2]}")
⋮----
# ============================================================
# Step 3a: Reconstruct clean slab (before adsorption)
⋮----
clean_sweeps = max(1, reconstruction_sweeps // 2)
⋮----
recon_clean_result = self.simulate_surface_reconstruction(
slab = recon_clean_result.lowest_energy_structure.atoms
⋮----
# Step 2: Adsorption (on reconstructed surface)
⋮----
adsorption_result = self.predict_adsorption_sites(
surface_with_adsorbate = adsorption_result.best_result.final_structure
⋮----
# Step 3b: Reconstruct with adsorbate present
⋮----
# Compute surface and adsorbate indices for proper constraint handling
n_slab = len(slab)
n_total = len(surface_with_adsorbate)
pre_mc_adsorbate_indices = list(range(n_slab, n_total))
slab_tags = slab.get_tags()
⋮----
pre_mc_surface_indices = [j for j in range(n_slab) if slab_tags[j] == 1]
⋮----
pre_mc_surface_indices = list(range(n_slab))
⋮----
z_coords = slab.positions[:, 2]
z_unique = sorted(set(round(z, 1) for z in z_coords), reverse=True)
top_2_z = set(z_unique[:2]) if len(z_unique) >= 2 else set(z_unique)
pre_mc_surface_indices = [
⋮----
reconstruction_result = self.simulate_surface_reconstruction(
# 使用重构后的最低能量结构（保留吸附物）
reconstructed_surface = reconstruction_result.lowest_energy_structure.atoms
# 识别吸附物原子索引，用于确定吸附位点
adsorbate_indices = [
⋮----
reconstruction_result = None
# 直接使用吸附后的结构
reconstructed_surface = surface_with_adsorbate
# 识别吸附物原子索引
adsorbate_indices = list(range(len(slab), len(reconstructed_surface)))
⋮----
# 步骤4: 反应路径分析（使用重构后的表面，通过替换吸附物进行）
pathway_result = self.analyze_reaction_pathway(
⋮----
num_sites=5,  # 重构后使用较少的位点数以节省时间
⋮----
current_adsorbate_indices=adsorbate_indices,  # 传递当前吸附物索引
⋮----
# 保存表面分析结果
analysis = SurfaceAnalysisResult(
⋮----
# 选择最优表面（基于能垒最低）
best_surface = None
⋮----
best_surface = min(
⋮----
# 步骤5: KMC模拟（可选）
kmc_result = None
⋮----
kmc_result = self.run_kmc_simulation(
⋮----
# 组装最终结果
result = GasSolidWorkflowResult(
⋮----
size: Tuple[int, int, int] = (2, 2, 1),  # 默认2x2 supercell
⋮----
"""
        从bulk结构生成slab

        参数:
            bulk: bulk结构
            miller: Miller指数
            layers: 层数
            vacuum: 真空层厚度 (Å)
            size: supercell大小 (a, b, c)，前两个是xy方向的扩展

        返回:
            slab结构
        """
⋮----
# 使用size参数生成supercell，确保表面足够大
slab = surface(bulk, miller, layers, vacuum=vacuum)
# 扩展xy方向
⋮----
slab = slab.repeat((size[0], size[1], 1))
⋮----
# 备用方案：返回简单的重复结构
slab = bulk.repeat((size[0] * 2, size[1] * 2, layers))
⋮----
# 便捷函数
⋮----
"""
    便捷函数：运行气固界面催化数字孪生

    参数:
        bulk_structure: bulk结构文件路径或ASE Atoms对象
        reaction_intermediates: 反应中间体列表（如["*O", "*OH", "*OOH", "*"]）
        initial_reactant: 初始反应物（如"*O"）
        temperature: 反应温度 (K)
        output_dir: 输出目录
        **kwargs: 其他参数传递给GasSolidDigitalTwin.run_complete_workflow

    返回:
        完整的工作流程结果
    """
# 配置日志
⋮----
# 检查是否有本地模型路径
fairchem_model_path = kwargs.pop('fairchem_model_path', None)
fairchem_model = kwargs.pop('fairchem_model', 'uma-s-1p1')
⋮----
# 如果没有提供本地路径，尝试使用项目中的模型
⋮----
project_model_path = os.path.join(
⋮----
fairchem_model_path = project_model_path
⋮----
# 初始化数字孪生
dt = GasSolidDigitalTwin(
⋮----
# 运行工作流程
result = dt.run_complete_workflow(
⋮----
# 示例：CO氧化反应
result = run_gas_solid_catalysis(
````

## File: camel_agents/mechanism_prompts.py
````python
"""Backward-compatibility re-exports.

All mechanism agent prompts now live in ``camel_agents.prompts`` alongside
Agent1-7 prompts.  This module re-exports for any code that imported from here.
"""
⋮----
from camel_agents.prompts import (  # noqa: F401
⋮----
__all__ = [
````

## File: camel_agents/mechanism_schemas.py
````python
"""CatDT universal mechanism search schema definitions.

These are CatDT's own representations for multi-pathway mechanism search,
independent of any particular backend (CARE, etc.).  They live between
Agent3 output and Agent4/5 input.
"""
⋮----
# ---------------------------------------------------------------------------
# Enums
⋮----
class ReactionOperationType(str, Enum)
⋮----
"""Types of elementary reaction operations."""
ADSORPTION = "adsorption"
DESORPTION = "desorption"
DISSOCIATION = "dissociation"
ASSOCIATION = "association"
HYDROGENATION = "hydrogenation"
DEHYDROGENATION = "dehydrogenation"
PCET = "pcet"                      # proton-coupled electron transfer
COUPLING = "coupling"              # C-C, C-N, etc.
REARRANGEMENT = "rearrangement"
ELEY_RIDEAL = "eley_rideal"
OTHER = "other"
⋮----
class StatePhase(str, Enum)
⋮----
"""Phase of a catalytic state species."""
GAS = "gas"
ADSORBED = "ads"
SURFACE = "surf"
SOLVATED = "solv"
⋮----
class EnergyBackendType(str, Enum)
⋮----
"""Backend used for quick free-energy estimation."""
THERMAL_UMA = "thermal_uma"        # UMA adsorption energy
ELECTRO_CHE = "electro_che"        # computational hydrogen electrode
HEURISTIC = "heuristic"            # rule-based estimate
EXTERNAL = "external"              # user-supplied value
⋮----
# State-level representations
⋮----
class CatalyticStateRecord(BaseModel)
⋮----
"""Lightweight, hashable record for one catalytic micro-state.

    This is intentionally richer than CARE's ``Intermediate``: it can encode
    surface identity, adsorption site type, coverage, and co-adsorbate
    configuration.
    """
state_id: str = Field(
species_label: str = Field(
phase: StatePhase = StatePhase.ADSORBED
elements: Dict[str, int] = Field(
surface_id: str = Field(
site_type: str = Field(
co_adsorbates: List[str] = Field(
coverage: float = Field(
structure_path: Optional[str] = Field(
extra: Dict[str, Any] = Field(
⋮----
class StateEnergyEstimate(BaseModel)
⋮----
"""Quick free-energy estimate for a single state."""
state_id: str
free_energy_eV: Optional[float] = None
adsorption_energy_eV: Optional[float] = None
backend: EnergyBackendType = EnergyBackendType.HEURISTIC
confidence: float = Field(default=0.0, ge=0.0, le=1.0)
details: Dict[str, Any] = Field(default_factory=dict)
⋮----
# Step-level representations
⋮----
class ElementaryStepCandidate(BaseModel)
⋮----
"""One candidate elementary reaction step between two catalytic states."""
step_id: str = Field(
reactant_state_id: str
product_state_id: str
operation_type: ReactionOperationType = ReactionOperationType.OTHER
bond_changes: str = Field(
estimated_barrier_eV: Optional[float] = None
estimated_reaction_energy_eV: Optional[float] = None
source_backend: str = Field(
⋮----
extra: Dict[str, Any] = Field(default_factory=dict)
⋮----
# Pathway-level representations
⋮----
class PathwayStateEntry(BaseModel)
⋮----
"""One state node in a candidate pathway."""
index: int
state: CatalyticStateRecord
energy_estimate: Optional[StateEnergyEstimate] = None
⋮----
class PathwayStepEntry(BaseModel)
⋮----
"""One step edge in a candidate pathway."""
⋮----
step: ElementaryStepCandidate
reactant_structure_path: Optional[str] = None
product_structure_path: Optional[str] = None
⋮----
class CandidatePathway(BaseModel)
⋮----
"""A single candidate reaction pathway (sequence of states + steps)."""
pathway_id: str
description: str = ""
states: List[PathwayStateEntry] = Field(default_factory=list)
steps: List[PathwayStepEntry] = Field(default_factory=list)
total_free_energy_change_eV: Optional[float] = None
max_step_energy_eV: Optional[float] = None
is_retained: bool = True
prune_reason: str = ""
export_dir: Optional[str] = None
⋮----
class MultiPathwaySearchResult(BaseModel)
⋮----
"""Top-level output of the mechanism search stage."""
context_summary: str = ""
all_pathways: List[CandidatePathway] = Field(default_factory=list)
retained_pathways: List[str] = Field(
pruned_pathways: List[str] = Field(
search_stats: Dict[str, Any] = Field(default_factory=dict)
export_manifest: Dict[str, Any] = Field(default_factory=dict)
⋮----
# Context / Config
⋮----
class MechanismContext(BaseModel)
⋮----
"""Structured context for the mechanism search stage.

    Created from Agent3 output + user-supplied reaction specification.
    """
bulk_formula: str = ""
surface_id: str = ""
surface_facet: str = ""
initial_state: str = Field(
target_state: str = Field(
environment: Dict[str, Any] = Field(
known_intermediates: List[str] = Field(
constraints: str = ""
surface_structure_path: Optional[str] = None
clean_slab_path: Optional[str] = None
surface_indices: List[int] = Field(default_factory=list)
adsorbate_indices: List[int] = Field(default_factory=list)
fixed_adsorption_site: Optional[List[float]] = Field(
⋮----
class ExplorationMode(str, Enum)
⋮----
"""Which exploration module(s) to run."""
AGENT_GUIDED = "agent_guided"      # LLM recommends pathway → fast evaluate
SYSTEMATIC = "systematic"          # CARE-style CRN build from scratch
BOTH_SEQUENTIAL = "both_sequential"  # agent-guided first, then systematic
BOTH_PARALLEL = "both_parallel"    # run both, merge results
⋮----
class MechanismSearchConfig(BaseModel)
⋮----
"""Tunable parameters for the search engine."""
exploration_mode: ExplorationMode = Field(
max_depth: int = Field(default=8, ge=1)
beam_width: int = Field(default=8, ge=1)
delta_keep_eV: float = Field(
delta_prune_eV: float = Field(
max_state_evaluations: int = Field(default=40, ge=1)
energy_backend: EnergyBackendType = EnergyBackendType.THERMAL_UMA
enable_care_backend: bool = True
enable_generic_ops_backend: bool = True
enable_llm_expansion: bool = False
⋮----
# Module exports
⋮----
__all__ = [
````

## File: camel_agents/mechanism_search.py
````python
"""Mechanism search stage runner.

Orchestrates Agent M1 and M2 (or runs tool-only mode) between Agent3 and
Agent4/5 in the CatDT workflow.
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
class MechanismStageRunner
⋮----
"""Runs the mechanism search stage after Agent3 completes.

    Can operate in two modes:
    1. **Agent mode**: Agent M1 and M2 coordinate via LLM + tool calls.
    2. **Tool-only mode**: Deterministic pipeline, no LLM involved.
       (Default for Phase 1 — simpler and more predictable.)
    """
⋮----
def __init__(self, tools: Any = None, use_agents: bool = False)
⋮----
"""Execute the mechanism search stage.

        Args:
            state: WorkflowState after Agent3 completion
            config: CatDTConfig with mechanism_* fields

        Returns:
            MultiPathwaySearchResult-compatible dict
        """
⋮----
# Extract context from workflow state
initial_state = ""
target_state = ""
known_intermediates: List[str] = []
⋮----
intermediates = list(state.reaction_context.intermediates or [])
⋮----
initial_state = intermediates[0]
target_state = intermediates[-1]
known_intermediates = intermediates[1:-1]
⋮----
surface_id = ""
surface_facet = ""
⋮----
# Try to extract facet from path, e.g. "slab_111.vasp"
m = re.search(r"(\d{3,4})", str(state.surface_path))
⋮----
surface_facet = m.group(1)
⋮----
# Step 1: Initialize context
⋮----
ctx_result = self.tools.initialize_mechanism_context(
⋮----
# Step 2: Generate candidate steps
candidates_result = self.tools.generate_candidate_steps(
⋮----
# Step 3: Run search
exploration_mode = getattr(config, "mechanism_exploration_mode", "agent_guided")
search_result = self.tools.run_mechanism_search(
⋮----
# Step 4: Export and prepare handoff
shortlist = self.tools.extract_pathway_shortlist(
⋮----
"""Convert mechanism search output to Agent4/5-ready payload.

        Returns:
            Dict with ``pathways`` list, each containing ``intermediates``
            and ``steps`` that Agent4/5 can consume.
        """
shortlist = search_output.get("shortlist", {})
payloads = shortlist.get("agent45_payloads", [])
⋮----
__all__ = ["MechanismStageRunner"]
````

## File: camel_agents/parallel_pathway.py
````python
"""
并行路径探索器 (ParallelPathwayExplorer)

支持：
1. 并行生成多个候选路径
2. 并行计算多个 NEB 步骤
3. 依赖分析和任务调度
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
@dataclass
class ParallelTask
⋮----
"""并行任务"""
task_id: str
task_type: str  # "pathway_design", "neb", "validation"
inputs: Dict[str, Any]
depends_on: List[str]  # 依赖的任务 ID
priority: int = 0
estimated_time: float = 0.0  # 预估时间（秒）
⋮----
@dataclass
class ParallelTaskResult
⋮----
"""并行任务结果"""
⋮----
success: bool
result: Any
error: Optional[str] = None
start_time: Optional[datetime] = None
end_time: Optional[datetime] = None
⋮----
@property
    def duration(self) -> float
⋮----
class DependencyGraph
⋮----
"""依赖图 - 管理任务间的依赖关系"""
⋮----
def __init__(self)
⋮----
self.dependencies: Dict[str, Set[str]] = defaultdict(set)  # task -> its dependencies
self.dependents: Dict[str, Set[str]] = defaultdict(set)    # task -> tasks that depend on it
⋮----
def add_task(self, task: ParallelTask)
⋮----
"""添加任务"""
⋮----
def get_ready_tasks(self, completed: Set[str]) -> List[ParallelTask]
⋮----
"""获取可以执行的任务（所有依赖已完成）"""
ready = []
⋮----
# 按优先级排序
⋮----
def get_parallel_groups(self) -> List[List[str]]
⋮----
"""
        将任务分组，每组内的任务可以并行执行
        
        Returns:
            任务 ID 的分组列表
        """
completed = set()
groups = []
remaining = set(self.tasks.keys())
⋮----
# 找到所有依赖已完成的任务
group = []
⋮----
# 有循环依赖
⋮----
def visualize(self) -> str
⋮----
"""生成依赖图的文本可视化"""
lines = ["Dependency Graph:"]
⋮----
deps = ", ".join(self.dependencies[task_id]) or "None"
⋮----
class ParallelPathwayExplorer
⋮----
"""
    并行路径探索器
    
    支持：
    1. 并行生成多个候选反应路径
    2. 并行计算 NEB（独立步骤）
    3. 智能任务调度
    """
⋮----
def __init__(self, max_workers: int = 4)
⋮----
"""
        并行探索多个反应路径
        
        Parameters:
            surface: 表面结构
            reaction_description: 反应描述
            pathway_design_fn: 路径设计函数 (surface, description) -> pathway
            n_pathways: 要生成的路径数量
            validation_fn: 验证函数 (pathway) -> validation_result
        
        Returns:
            多个路径设计结果，按 confidence 排序
        """
⋮----
# 提交并行任务
pathways = []
⋮----
# 提交多个路径设计任务
futures = []
⋮----
future = executor.submit(
⋮----
seed=i  # 不同的随机种子产生不同路径
⋮----
# 收集结果
⋮----
pathway = future.result(timeout=300)  # 5分钟超时
⋮----
# 验证（如果提供了验证函数）
⋮----
validated = []
⋮----
futures = {
⋮----
pathway = futures[future]
⋮----
validation = future.result()
⋮----
pathways = validated
⋮----
# 按 confidence 排序
⋮----
"""使用特定随机种子设计路径"""
⋮----
"""
        并行计算多个 NEB 步骤
        
        自动分析步骤依赖关系，并行计算独立的步骤。
        
        Parameters:
            steps: 步骤列表，每个包含 name, reactant, product
            neb_fn: NEB 计算函数 (step) -> neb_result
            max_workers: 并行工作线程数
        
        Returns:
            步骤名称 -> NEB 结果 的字典
        """
workers = max_workers or self.max_workers
⋮----
# 分析依赖关系
# 步骤 i 的结果可能作为步骤 i+1 的输入
# 但在这里我们假设所有步骤是独立的（可以同时计算）
⋮----
results = {}
completed = 0
failed = 0
⋮----
# 提交所有任务
future_to_step = {
⋮----
step = future_to_step[future]
step_name = step.get('name', 'unknown')
⋮----
result = future.result()
⋮----
"""
        使用依赖图执行任务
        
        Parameters:
            tasks: 任务列表
            execute_fn: 执行函数 (task) -> result
        
        Returns:
            任务 ID -> 结果 的字典
        """
# 构建依赖图
graph = DependencyGraph()
⋮----
# 获取并行组
groups = graph.get_parallel_groups()
⋮----
# 执行
⋮----
group_tasks = [graph.tasks[tid] for tid in group]
⋮----
future_to_task = {
⋮----
task = future_to_task[future]
result = ParallelTaskResult(
⋮----
"""
        从多个路径中选择最佳路径
        
        Parameters:
            pathways: 路径列表
            criteria: 选择标准 ("confidence", "lowest_barrier", "fewest_steps")
        
        Returns:
            最佳路径
        """
⋮----
# 需要 NEB 结果
def max_barrier(p)
⋮----
barriers = [
⋮----
# =============================================================================
# 便捷函数
⋮----
"""
    并行计算多个反应步骤的能垒
    
    Example:
        >>> steps = [
        ...     {"name": "CO_to_CHO", "reactant": atoms1, "product": atoms2},
        ...     {"name": "CHO_to_CH2O", "reactant": atoms2, "product": atoms3},
        ... ]
        >>> barriers = parallel_compute_barriers(steps, neb_calculator)
        >>> print(barriers)
        {"CO_to_CHO": 0.84, "CHO_to_CH2O": 1.23}
    """
explorer = ParallelPathwayExplorer(max_workers=max_workers)
⋮----
def neb_fn(step)
⋮----
"""
    比较多个路径的性能指标
    
    Returns:
        比较报告
    """
report = {
⋮----
values = [p.get('confidence', 0) for p in pathways]
⋮----
values = [len(p.get('steps', [])) for p in pathways]
⋮----
# 综合排名
scores = []
⋮----
score = p.get('confidence', 0) * 0.5  # confidence 权重 50%
score += (1.0 / max(len(p.get('steps', [])), 1)) * 0.3  # 步数少加分
````

## File: camel_agents/pathway_utils.py
````python
logger = logging.getLogger(__name__)
⋮----
DEFAULT_RENDERER = "tachyon"
DEFAULT_CLEARANCE = 1.1
⋮----
def get_default_her_description() -> str
⋮----
def get_default_nrr_description() -> str
⋮----
def align_added_atoms_to_anchor(positions: List[List[float]], anchor_position: List[float]) -> List[List[float]]
⋮----
anchor = [float(anchor_position[0]), float(anchor_position[1]), float(anchor_position[2])]
base = positions[0]
aligned: List[List[float]] = []
⋮----
dx = float(pos[0]) - float(base[0])
dy = float(pos[1]) - float(base[1])
dz = float(pos[2]) - float(base[2])
⋮----
# First atom sits exactly above anchor in x/y, retains original z
⋮----
surface_indices = list(range(len(slab)))
⋮----
top_idx = max(surface_indices, key=lambda idx: slab.positions[idx][2])
⋮----
min_allowed = float(surface_max_z) + float(clearance)
adjusted: List[List[float]] = []
⋮----
z = min_allowed
⋮----
resolved: List[int] = []
⋮----
idx_raw = raw.get("index") if isinstance(raw, dict) else raw
⋮----
idx = int(idx_raw)
⋮----
product = structure.copy()
adsorbate_indices = list(product.info.get("adsorbate_indices", []))
⋮----
# Add atoms
⋮----
symbol = (add_spec.species or "").strip()[:2]
⋮----
symbol = symbol[0].upper() + (symbol[1].lower() if len(symbol) > 1 else "")
⋮----
position = [float(add_spec.position[0]), float(add_spec.position[1]), float(add_spec.position[2])]
⋮----
position = [0.0, 0.0, float(product.positions[:, 2].max() + 1.5)]
⋮----
# Remove atoms (protect surface-only structures in local-index mode)
to_remove = _resolve_remove_indices(step.atoms_to_remove, adsorbate_indices, use_adsorbate_local_indices)
⋮----
to_remove = []
⋮----
remove_set = set(i for i in to_remove if 0 <= i < len(product))
keep_mask = [i not in remove_set for i in range(len(product))]
index_map: Dict[int, int] = {}
next_idx = 0
⋮----
product = product[keep_mask]
adsorbate_indices = [index_map[i] for i in adsorbate_indices if i in index_map]
⋮----
ads_set = set(adsorbate_indices)
⋮----
labels = list(intermediates)
energy_list = [float(e) for e in energies]
barrier_map: Dict[str, float] = {}
⋮----
key = f"{intermediates[i]} -> {intermediates[i + 1]}"
⋮----
@dataclass
class PathwayRunState
⋮----
reaction_description: str
anchor_position: List[float]
surface_indices: List[int]
adsorbate_indices: List[int]
⋮----
class PathwayAgentsRunner
⋮----
"""Compatibility runner placeholder.

    Full pathway orchestration is now handled by `CatDTCamelWorkflow`.
    """
⋮----
def __init__(self, *args: Any, **kwargs: Any)
⋮----
_ = args, kwargs
⋮----
def run(self, state: PathwayRunState, base_structure: Atoms) -> PathwayRunState
⋮----
_ = base_structure
````

## File: camel_agents/policy.py
````python
"""Evolvability policy for CAMEL workflow strategy adaptation."""
⋮----
logger = logging.getLogger(__name__)
⋮----
class EvolvablePolicy
⋮----
"""Simple bandit-style evolvability policy updated from workflow rewards."""
⋮----
def _load(self) -> None
⋮----
payload = json.loads(self.policy_path.read_text(encoding="utf-8"))
q_values = payload.get("q_values", {})
visits = payload.get("visits", {})
⋮----
rng_state = str(payload.get("rng_state", "") or "").strip()
⋮----
def _save(self) -> None
⋮----
payload = {
⋮----
def choose_strategy(self) -> str
⋮----
def update(self, strategy: str, reward: float) -> None
⋮----
old = self.q_values[strategy]
````

## File: camel_agents/prompts_zh_backup.py
````python
"""
CatDT Agent Prompts (CAMEL-ready)

设计原则：
1. 提示词完全通用，不包含针对某一反应的规则或示例。
2. Agent4/Agent5 通过多轮反馈协作（最多 N 轮）逐步修正。
3. 不预设固定元素增删规则，所有增删改由 LLM 基于输入自行判断。
4. 完整结构坐标仅在首次输入一次，后续只给吸附分子坐标。
"""
⋮----
@dataclass
class AgentRole
⋮----
name: str
role: str
goal: str
backstory: str
⋮----
def to_dict(self) -> Dict[str, str]
⋮----
class TaskPrompts
⋮----
"""任务提示词模板"""
⋮----
@staticmethod
    def reaction_context_parsing(reaction_description: str) -> str
⋮----
"""Agent4 路径端点修正提示词。

        Args:
            reaction_context: 反应上下文 dict (reaction_type, intermediates, constraints)
            surface_info: 表面概要信息（表面顶部 z 坐标、晶胞参数等）
            baseline_steps: 每步吸附分子坐标与元素差异（compact 格式）
            staging_plan_text: 程序预计算的暂驻方案（元素差异 + 建议坐标）
            feedback: 上一轮 Agent5 反馈
            history: 迭代历史摘要
        """
feedback_section = f"\n【上一轮反馈（含前一轮吸附分子坐标）】\n{feedback}\n" if feedback else ""
history_section = f"\n【迭代历史】\n{history}\n" if history else ""
staging_section = staging_plan_text if staging_plan_text else "(无需暂驻)"
memento_section = memento_context if memento_context else "(no retrieved memento cases)"
knowledge_section = knowledge_context if knowledge_context else "(no retrieved knowledge items)"
skill_section = skill_context if skill_context else "(no retrieved skill items)"
⋮----
"""Agent5 端点验证提示词。

        Args:
            steps_details: 每步摘要（step/formula/adds/removes）
            structures_payload: 每步吸附分子坐标（compact 格式）
            precheck_summary: 程序预检概览
            precheck_issues: 程序预检问题明细（fatal/warning）
            staging_plan_text: 预计算暂驻方案参考
        """
precheck_text = precheck_summary if precheck_summary else "(none)"
precheck_issue_text = precheck_issues if precheck_issues else "(none)"
⋮----
# ----- Agent 1-3, 6-7 简短提示词 -----
⋮----
@staticmethod
    def agent6_run_neb(step_structures: List[Dict[str, Any]]) -> str
⋮----
step_names = [s.get("name", "unknown") for s in step_structures]
⋮----
# ----- Agent Role Definitions -----
⋮----
AGENT1_STRUCTURE_INITIALIZER = AgentRole(
⋮----
AGENT2_ADSORPTION_PREDICTOR = AgentRole(
⋮----
AGENT3_RECONSTRUCTION_SIMULATOR = AgentRole(
⋮----
AGENT4_PATHWAY_DESIGNER = AgentRole(
⋮----
AGENT5_PATHWAY_VALIDATOR = AgentRole(
⋮----
AGENT6_NEB_RUNNER = AgentRole(
⋮----
AGENT7_REPORT_GENERATOR = AgentRole(
⋮----
def get_agent_role(agent_name: str) -> Optional[AgentRole]
⋮----
roles = {
⋮----
def get_all_agent_roles() -> Dict[str, AgentRole]
````

## File: camel_agents/prompts.py
````python
"""
CatDT Agent Prompts (CAMEL-ready)

Design principles:
1. Prompts are fully generic — no reaction-specific rules or examples.
2. Agent4/Agent5 collaborate via multi-round feedback (up to N rounds).
3. No hardcoded element add/remove rules — LLM decides based on input.
4. Full structure coordinates given once; subsequent rounds only show adsorbate coords.

Backup of previous Chinese version: prompts_zh_backup.py
"""
⋮----
@dataclass
class AgentRole
⋮----
name: str
role: str
goal: str
backstory: str
⋮----
def to_dict(self) -> Dict[str, str]
⋮----
class TaskPrompts
⋮----
"""任务提示词模板"""
⋮----
@staticmethod
    def reaction_context_parsing(reaction_description: str) -> str
⋮----
"""Agent4 路径端点修正提示词。

        Args:
            reaction_context: 反应上下文 dict (reaction_type, intermediates, constraints)
            surface_info: 表面概要信息（表面顶部 z 坐标、晶胞参数等）
            baseline_steps: 每步吸附分子坐标与元素差异（compact 格式）
            feedback: 上一轮 Agent5 反馈
            history: 迭代历史摘要
        """
feedback_section = f"\n【上一轮反馈（含前一轮吸附分子坐标）】\n{feedback}\n" if feedback else ""
history_section = f"\n【迭代历史】\n{history}\n" if history else ""
memento_section = memento_context if memento_context else "(no retrieved memento cases)"
knowledge_section = knowledge_context if knowledge_context else "(no retrieved knowledge items)"
skill_section = skill_context if skill_context else "(no retrieved skill items)"
⋮----
"""Agent5 端点验证提示词。

        Args:
            steps_details: 每步摘要（step/formula/adds/removes）
            structures_payload: 每步吸附分子坐标（compact 格式）
            precheck_summary: 程序预检概览
            precheck_issues: 程序预检问题明细（fatal/warning）
        """
precheck_text = precheck_summary if precheck_summary else "(none)"
precheck_issue_text = precheck_issues if precheck_issues else "(none)"
⋮----
# ----- Mechanism search prompts -----
⋮----
"""Prompt LLM to recommend multiple competing reaction pathways."""
constraint_line = f"\nConstraints: {constraints}" if constraints else ""
⋮----
"""Agent M1 task prompt."""
⋮----
"""Agent M2 task prompt."""
⋮----
# ----- Agent 1-3, 6-7 简短提示词 -----
⋮----
@staticmethod
    def agent6_run_neb(step_structures: List[Dict[str, Any]]) -> str
⋮----
step_names = [s.get("name", "unknown") for s in step_structures]
⋮----
# ----- Agent Role Definitions -----
⋮----
AGENT1_STRUCTURE_INITIALIZER = AgentRole(
⋮----
AGENT2_ADSORPTION_PREDICTOR = AgentRole(
⋮----
AGENT3_RECONSTRUCTION_SIMULATOR = AgentRole(
⋮----
AGENT4_PATHWAY_DESIGNER = AgentRole(
⋮----
AGENT5_PATHWAY_VALIDATOR = AgentRole(
⋮----
AGENT6_NEB_RUNNER = AgentRole(
⋮----
AGENT7_REPORT_GENERATOR = AgentRole(
⋮----
# --- Mechanism search agents (between Agent3 and Agent4/5) ---
⋮----
MECHANISM_CONTEXT_ROUTING_AGENT = AgentRole(
⋮----
MECHANISM_SEARCH_TRIAGE_AGENT = AgentRole(
⋮----
def get_agent_role(agent_name: str) -> Optional[AgentRole]
⋮----
roles = {
⋮----
def get_all_agent_roles() -> Dict[str, AgentRole]
````

## File: camel_agents/runtime.py
````python
"""CAMEL runtime primitives for agent task execution."""
⋮----
@dataclass
class TaskOutput
⋮----
raw: str = ""
json_dict: Dict[str, Any] = field(default_factory=dict)
pydantic: Optional[Any] = None
tool_call_record: Optional[Any] = None
⋮----
@dataclass
class Task
⋮----
description: str
expected_output: str
agent: Any
output_pydantic: Optional[Any] = None
result_handler: Optional[Any] = None
allow_tool_calls: Optional[bool] = None
⋮----
class CamelWorkflowAgent
⋮----
@staticmethod
    def _prepare_tools(tools: Optional[List[FunctionTool]]) -> List[FunctionTool]
⋮----
"""Clone tools and freeze validated schemas to avoid runtime schema drift."""
prepared: List[FunctionTool] = []
⋮----
tool_name = tool.get_function_name()
⋮----
schema = copy.deepcopy(tool.get_openai_tool_schema())
cloned_tool = FunctionTool(tool.func, openai_tool_schema=schema)
cloned_tool.get_openai_tool_schema = (  # type: ignore[method-assign]
⋮----
cloned_tool.get_openai_function_schema = (  # type: ignore[method-assign]
⋮----
system_message = (
prepared_tools = self._prepare_tools(tools or [])
# No max_iteration — let the loop run until LLM stops calling tools.
# Timeout is controlled by step_timeout (set per-agent).
⋮----
@staticmethod
    def _extract_json(text: str) -> Dict[str, Any]
⋮----
stripped = text.strip()
block = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", stripped, re.DOTALL | re.IGNORECASE)
candidate = block.group(1) if block else stripped
⋮----
parsed = json.loads(candidate)
⋮----
# Try full-text direct parse first
⋮----
parsed = json.loads(stripped)
⋮----
# Fallback: scan all balanced {...} segments and parse from tail.
segments: List[str] = []
start = None
depth = 0
in_string = False
escaped = False
⋮----
escaped = True
⋮----
in_string = True
⋮----
start = idx
⋮----
parsed = json.loads(seg)
⋮----
@staticmethod
    def _model_validate(output_model: Any, payload: Dict[str, Any]) -> Any
⋮----
@staticmethod
    def _is_length_finish_error(exc: Exception) -> bool
⋮----
text = str(exc or "")
exc_name = exc.__class__.__name__
⋮----
@staticmethod
    def _is_response_format_parse_error(exc: Exception) -> bool
⋮----
lowered = text.lower()
⋮----
# Root-cause fix:
# OpenAI-compatible structured parse (`response_format`) is unstable across SDK/model variants
# in our runtime stack (e.g., typing_extensions.Required / AnnotatedAlias transform failures).
# Keep a single generic pathway: plain text generation + local JSON extraction + local pydantic validation.
# Can be re-enabled for experiments by setting CATDT_ENABLE_SERVER_STRUCTURED_PARSE=1.
enable_server_structured_parse = os.getenv("CATDT_ENABLE_SERVER_STRUCTURED_PARSE", "0").strip() in {
response_format = (
before_sig = getattr(self.chat_agent, "_last_tool_call_signature", None)
⋮----
def _step_with_optional_tools(rf: Optional[Any])
⋮----
tool_dict = getattr(self.chat_agent, "tool_dict", {}) or {}
removed_tools = list(tool_dict.values())
removed_names = [tool.get_function_name() for tool in removed_tools]
⋮----
response = _step_with_optional_tools(response_format)
⋮----
# Structured parsing can fail when model output is truncated at max tokens.
# Retry once without server-side response_format, then parse JSON locally.
⋮----
response = _step_with_optional_tools(None)
⋮----
after_sig = getattr(self.chat_agent, "_last_tool_call_signature", None)
tool_call_record = None
⋮----
tool_call_record = getattr(self.chat_agent, "_last_tool_call_record", None)
content = ""
⋮----
content = getattr(response.msgs[-1], "content", "") or ""
⋮----
# If tool calls were made but the final response has no text content,
# do a follow-up call WITHOUT tools to force the LLM to output JSON.
⋮----
follow_up = (
⋮----
follow_up_response = self.chat_agent.step(follow_up, response_format=None)
⋮----
content = getattr(follow_up_response.msgs[-1], "content", "") or ""
⋮----
payload = self._extract_json(content)
parsed = self._model_validate(output_model, payload)
````

## File: camel_agents/schemas.py
````python
"""Shared schema definitions for CAMEL CatDT workflow."""
⋮----
class ReactionContext(BaseModel)
⋮----
reaction_type: str
initial_adsorbate: str
intermediates: List[str]
constraints: str = ""
⋮----
class AtomAddSpec(BaseModel)
⋮----
species: str = ""
element: str = ""
symbol: str = ""
position: Optional[List[float]] = None
reactant_position: Optional[List[float]] = None
product_position: Optional[List[float]] = None
reason: str = ""
⋮----
class AtomModifySpec(BaseModel)
⋮----
index: int
⋮----
class PathwayStepSpec(BaseModel)
⋮----
step_name: str
reactant_formula: str
product_formula: str
⋮----
atoms_to_add: List[AtomAddSpec] = Field(default_factory=list)
atoms_to_modify: List[AtomModifySpec] = Field(default_factory=list)
atoms_to_remove: List[Any] = Field(default_factory=list)
confidence: float = 0.0
⋮----
class PathwayDesign(BaseModel)
⋮----
pathway_name: str
description: str
overall_reaction: str
steps: List[PathwayStepSpec]
⋮----
class ValidationReport(BaseModel)
⋮----
status: str
issues: List[str] = Field(default_factory=list)
feedback: str = ""
⋮----
@dataclass
class WorkflowState
⋮----
run_id: str = ""
output_base_dir: str = ""
reaction_description: str = ""
surface_indices: Optional[List[int]] = None
adsorbate_indices: Optional[List[int]] = None
surface_atom_count: Optional[int] = None
reaction_context: Optional[ReactionContext] = None
surface_path: Optional[str] = None
clean_slab_path: Optional[str] = None
surface_with_adsorbate_path: Optional[str] = None
reconstructed_surface_path: Optional[str] = None
intermediates: List[str] = field(default_factory=list)
pathway_design: Optional[PathwayDesign] = None
validation_report: Optional[ValidationReport] = None
step_structures: List[Dict[str, Any]] = field(default_factory=list)
tool_baseline_steps: List[Dict[str, Any]] = field(default_factory=list)
memento_retrieval: Dict[str, Any] = field(default_factory=dict)
knowledge_retrieval: Dict[str, Any] = field(default_factory=dict)
skill_retrieval: Dict[str, Any] = field(default_factory=dict)
memento_query_text: str = ""
energy_gate_report: Dict[str, Any] = field(default_factory=dict)
adsorption_energies: Dict[str, float] = field(default_factory=dict)
adsorbate_energies: Dict[str, float] = field(default_factory=dict)
neb_results: Dict[str, Any] = field(default_factory=dict)
last_feedback: str = ""
iteration_history: List[Dict[str, Any]] = field(default_factory=list)
⋮----
# Multi-pathway mechanism search fields (populated when enable_mechanism_search=True)
mechanism_context: Optional[Dict[str, Any]] = None
candidate_pathways: List[Dict[str, Any]] = field(default_factory=list)
retained_pathways: List[Dict[str, Any]] = field(default_factory=list)
pruned_pathways: List[Dict[str, Any]] = field(default_factory=list)
pathway_manifests: Dict[str, Any] = field(default_factory=dict)
mechanism_search_result: Optional[Dict[str, Any]] = None
⋮----
# Multi-facet fields (populated when multi_facet=True)
all_surface_paths: Dict[str, str] = field(default_factory=dict)
all_area_fractions: Dict[str, float] = field(default_factory=dict)
all_miller_indices: Dict[str, tuple] = field(default_factory=dict)
facet_results: List[Any] = field(default_factory=list)  # List[FacetResult]
aggregated_tof: Optional[float] = None
aggregated_production_rates: Dict[str, float] = field(default_factory=dict)
current_facet_id: Optional[str] = None
⋮----
@dataclass
class FacetResult
⋮----
"""Per-facet pipeline output for multi-facet Wulff-weighted aggregation."""
facet_id: str = ""
miller_index: Optional[tuple] = None
area_fraction: float = 0.0
⋮----
pathway_result_pkl: Optional[str] = None
kmc_result_pkl: Optional[str] = None
tof: float = 0.0
production_rates: Dict[str, float] = field(default_factory=dict)
neb_summary: Dict[str, Any] = field(default_factory=dict)
error: Optional[str] = None
⋮----
class CatDTConfig(BaseModel)
⋮----
run_id: str = Field(
output_base_dir: str = Field(default="output/catdt_workflow")
⋮----
bulk_structure_path: Optional[str] = None
initial_surface_path: Optional[str] = None
reaction_description: Optional[str] = None
⋮----
top_n_surfaces: int = Field(3, ge=1)
num_adsorption_sites: int = Field(5, ge=1)
enable_llm_adsorption_review: bool = False
llm_adsorption_review_context: Optional[str] = None
⋮----
adsorbate_elements_for_mc: Optional[List[str]] = None
mc_temperature_k: Optional[float] = Field(None, gt=0)
mc_total_sweeps: int = Field(100, ge=1)
clean_slab_path_for_mc: Optional[str] = None
⋮----
# Two-step reconstruction: Agent3a(clean) → Agent2 → Agent3b(with adsorbate)
mc_two_step_reconstruction: bool = Field(
mc_clean_sweeps: Optional[int] = Field(
⋮----
# Energy model for VSSR-MC surface reconstruction
mc_energy_model: str = Field(
⋮----
# Electrochemical parameters (liquid-solid interface)
potential_she: Optional[float] = Field(
ph: Optional[float] = Field(
⋮----
# Materials Project API key (for Pourbaix diagram data)
mp_api_key: Optional[str] = Field(
⋮----
# Multi-facet Wulff-weighted microkinetics
multi_facet: bool = Field(
multi_facet_top_n: int = Field(
⋮----
calculate_barriers: bool = True
num_pathway_sites: int = Field(3, ge=1)
neb_n_frames: int = Field(5, ge=3)
neb_max_steps: int = Field(80, ge=1)
adsorption_relax_max_steps: int = Field(120, ge=1)
⋮----
kmc_temperature_k: Optional[float] = Field(None, gt=0)
gas_pressures: Optional[Dict[str, float]] = None
⋮----
# Mechanism search stage (Agent3 → Agent4/5 upstream)
enable_mechanism_search: bool = Field(
mechanism_exploration_mode: str = Field(
mechanism_max_depth: int = Field(8, ge=1)
mechanism_beam_width: int = Field(8, ge=1)
mechanism_delta_keep: float = Field(
mechanism_delta_prune: float = Field(
mechanism_max_state_evaluations: int = Field(40, ge=1)
⋮----
class Config
⋮----
arbitrary_types_allowed = True
⋮----
@classmethod
    def create_from_args(cls, **kwargs)
⋮----
@classmethod
    def from_kwargs(cls, **kwargs)
````

## File: camel_agents/tools.py
````python
"""CAMEL tool registry facade for CatDT."""
⋮----
class CatDTTools(
⋮----
"""Unified CatDT tool facade; implementations are split across tooling modules."""
⋮----
def _init_workflow_bridge(self) -> None
⋮----
"""Expose workflow helper methods via tools (for orchestration delegation)."""
⋮----
class _WorkflowToolBridge(WorkflowContextMixin, WorkflowGeometryMixin)
⋮----
bridge = _WorkflowToolBridge()
⋮----
def __getattr__(self, name: str) -> Any
⋮----
bridge = self.__dict__.get("_workflow_bridge")
⋮----
def build_steps_payload_tool(self, step_structures: list[dict] | None = None) -> str
⋮----
"""LLM-facing wrapper for step payload generation."""
⋮----
def programmatic_validation_tool(self, step_structures: list[dict] | None = None) -> dict
⋮----
"""LLM-facing wrapper for generic geometric validation."""
⋮----
def get_step_element_deltas_tool(self) -> list[dict]
⋮----
"""Query per-step element deltas: what atoms must be added to which side for NEB.

        Returns a list of steps, each with: step_name, to_add_to_reactant,
        to_add_to_product, total_staged.  Call this FIRST to plan your work.
        """
⋮----
"""Suggest positions for ALL staged atoms across ALL steps in one call.

        Pass a list of requests, each: {"step_name": "...", "element": "H", "side": "reactant"}.
        Returns suggested [x,y,z] positions near chemically relevant anchor atoms.
        Only call this if you need guidance — you can also propose positions yourself.
        """
⋮----
"""Validate ALL proposed atom positions across ALL steps in one call.

        Pass a list of steps, each: {"step_name": "...", "atoms": [{"species": "H", "side": "reactant", "position": [x,y,z]}]}.
        Returns per-step and per-atom diagnostics: distance checks, overlap checks, and suggestions.
        """
⋮----
"""UMA endpoint audit + constrained staged-atom relaxation for Agent4/5 loop."""
⋮----
sample = step_structures[0]
⋮----
def get_agent45_parametric_retriever_status_tool(self) -> dict
⋮----
# ---- Mechanism search tool wrappers ----
⋮----
"""Initialize mechanism search context from reaction specification.

        Args:
            fixed_adsorption_site: Optional [x,y,z] site from Agent2/3.
                If set, all intermediates are evaluated at this fixed site,
                skipping the expensive adsorption site search.
        """
⋮----
"""Generate candidate elementary steps from enabled backends."""
⋮----
"""Run mechanism search using selected exploration module(s).

        Args:
            exploration_mode: 'agent_guided', 'systematic', 'both_sequential', 'both_parallel'
        """
⋮----
"""Export retained pathways and prepare Agent4/5 handoff payload."""
⋮----
def to_camel_tools(self) -> Dict[str, FunctionTool]
⋮----
"""Expose CatDT tool methods as CAMEL `FunctionTool` registry."""
⋮----
# Mechanism search tools
⋮----
__all__ = [
````

## File: camel_agents/visualized_digital_twin.py
````python
"""
Gas-Solid Digital Twin with Integrated Visualization - COMPLETE VERSION

完整实现所有6类可视化：
1. ✅ SurFF 表面生成后 → 结构图片
2. ✅ AdsorbDiff 吸附位点后 → 结构图片
3. ✅ MC 表面重构 → 轨迹 GIF
4. ✅ 反应路径 → 完整过程 GIF
5. ✅ 能垒计算 → 自由能图
6. ✅ HTML摘要

修复内容：
- 完整实现02_adsorption可视化
- 完整实现03_mc_reconstruction可视化
- 完整实现05_energy_diagram可视化
- 修复反应路径顺序问题
- 添加自动扩胞逻辑

Author: Claude
Date: 2026-01-22
"""
⋮----
# Import visualization manager
⋮----
# Import existing digital twin
⋮----
logger = logging.getLogger(__name__)
⋮----
def check_adsorbate_distance(atoms: Atoms, threshold: float = 3.0) -> tuple
⋮----
"""
    检查吸附分子之间的距离是否过近

    Returns
    -------
    (is_too_close, min_distance, needs_expansion)
    """
⋮----
# 找出吸附分子（非Pt原子）
surface_atoms = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == 'Pt']
ads_atoms = [i for i in range(len(atoms)) if i not in surface_atoms]
⋮----
ads_pos = atoms.get_positions()[ads_atoms]
min_dist = pdist(ads_pos).min()
⋮----
is_too_close = min_dist < threshold
needs_expansion = is_too_close
⋮----
def expand_supercell(atoms: Atoms, factor: tuple = (2, 2, 1)) -> Atoms
⋮----
"""
    扩展supercell以避免吸附分子过近

    Parameters
    ----------
    atoms : Atoms
        原始结构
    factor : tuple
        扩展因子 (nx, ny, nz)

    Returns
    -------
    Atoms
        扩展后的结构
    """
⋮----
P = np.diag(factor)
expanded = make_supercell(atoms, P)
⋮----
class VisualizedGasSolidDigitalTwin(GasSolidDigitalTwin)
⋮----
"""
    带完整可视化的气固催化数字孪生系统

    扩展原有的 GasSolidDigitalTwin，在关键步骤自动生成所有可视化
    """
⋮----
# VSSR-MC model and electrochemical parameters (forwarded to parent)
⋮----
# 可视化参数
⋮----
# 扩胞参数
⋮----
# LLM agent审查参数
⋮----
# 可视化设置
⋮----
# 如果没有指定renderer，使用系统默认（按优先级：tachyon > povray > matplotlib）
⋮----
viz_renderer = 'tachyon'
⋮----
self.viz_manager = None  # 延迟初始化
⋮----
# 扩胞设置
⋮----
# LLM agent审查设置
⋮----
api_key = llm_api_key or os.getenv("OPENAI_API_KEY")
base_url = llm_base_url or os.getenv("OPENAI_BASE_URL")
⋮----
def _init_viz_manager(self, output_dir: str)
⋮----
"""初始化可视化管理器"""
⋮----
"""使用 OpenAI-compatible LLM 审查吸附结构合理性。"""
⋮----
chemical_formula = structure.get_chemical_formula()
positions = structure.get_positions()
symbols = structure.get_chemical_symbols()
cell = structure.cell.array
⋮----
prompt = f"""You are an expert in computational catalysis. Review this adsorbed structure and return JSON only.
⋮----
response = self.llm_client.chat.completions.create(
content = response.choices[0].message.content or ""
⋮----
json_match = re.search(r'\{.*\}', content, re.DOTALL)
⋮----
review_result = json.loads(json_match.group())
⋮----
review_result = {
⋮----
"""
        运行完整工作流并生成所有可视化

        Returns
        -------
        Dict
            包含计算结果和所有可视化文件路径的字典
        """
⋮----
# 初始化可视化管理器
⋮----
all_visualizations = {
⋮----
'kmc_dynamics': '',  # KMC表面动力学轨迹
⋮----
# 先调用父类的完整工作流（这会执行所有步骤包括KMC）
⋮----
workflow_result = super().run_complete_workflow(
⋮----
# 现在从结果中生成完整的可视化
⋮----
# 1. 可视化所有生成的表面 [01_surfaces]
⋮----
surface_files = {}
surface_energies = {}
⋮----
miller_str = ''.join(map(str, analysis.miller_index))
slab_file = os.path.join(output_dir, f"slab_{miller_str}.vasp")
⋮----
viz_paths = self.viz_manager.visualize_generated_surfaces(
⋮----
# 2. 可视化吸附位点 [02_adsorption]
⋮----
best = workflow_result.best_surface
miller_str = ''.join(map(str, best.miller_index))
⋮----
# 查找吸附位点文件
ads_dir = os.path.join(output_dir, f"surface_{miller_str}", "adsorption")
best_config_file = os.path.join(ads_dir, "best_config.vasp")
⋮----
# 创建一个字典，包含最佳吸附构型
ads_structures = {
⋮----
# 也可以加入diff轨迹的最终状态
traj_dir = os.path.join(ads_dir, "trajectories")
⋮----
for i in range(min(num_adsorption_sites, 5)):  # 最多5个
traj_file = os.path.join(traj_dir, f"diff_{i}", "0.traj")
⋮----
traj = Trajectory(traj_file)
⋮----
# 保存最终构型为vasp
final_config = os.path.join(ads_dir, f"config_{i:02d}.vasp")
⋮----
ads_energies = {k: 0.0 for k in ads_structures.keys()}  # TODO: get real energies
viz_paths = self.viz_manager.visualize_adsorption_sites(
⋮----
# 3. 可视化MC重构轨迹 [03_mc_reconstruction]
⋮----
# 从MC的CIF文件创建轨迹
recon_dir = os.path.join(output_dir, f"surface_{miller_str}", "reconstruction")
⋮----
cif_files = glob.glob(os.path.join(recon_dir, "inb_unrelaxed_slab_*.cif"))
⋮----
# 按照sweep顺序排序
⋮----
# 创建trajectory文件
traj_file = os.path.join(recon_dir, "trajectory.traj")
⋮----
traj = AseTrajectory(traj_file, 'w')
⋮----
for cif_file in cif_files[:20]:  # 最多20帧
⋮----
atoms = read(cif_file)
⋮----
mc_gif = self.viz_manager.visualize_mc_reconstruction(
⋮----
# 4. 可视化反应路径 [04_reaction_pathway]
⋮----
# 收集反应路径的结构文件（按正确的顺序）
pathway_dir = os.path.join(output_dir, f"surface_{miller_str}", "pathway", "adsorption")
⋮----
pathway_structures = {}
pathway_files = {}
⋮----
# 按照reaction_intermediates的顺序收集结构
⋮----
# 文件名是 "*CO_relaxed.vasp" 格式（星号是文件名的一部分）
inter_file = os.path.join(pathway_dir, f"{inter}_relaxed.vasp")
⋮----
atoms = read(inter_file)
⋮----
# 检查是否需要扩胞
final_structure = atoms
⋮----
expanded = expand_supercell(atoms, (2, 2, 1))
expanded_file = inter_file.replace('.vasp', '_2x2.vasp')
⋮----
final_structure = expanded
inter_file = expanded_file
⋮----
# LLM审查结构（如果启用）
⋮----
# 先生成单独的图片用于审查
temp_img_dir = os.path.join(output_dir, f"surface_{miller_str}", "pathway", "review_imgs")
⋮----
temp_img_path = os.path.join(temp_img_dir, f"{inter}_preview.png")
⋮----
# 快速生成可视化图片
⋮----
# 如果ASE write失败，使用viz_manager
⋮----
temp_img_path = self.viz_manager._render_structure(
⋮----
# 调用LLM审查
review_result = self._review_structure_with_llm(
⋮----
# 如果LLM发现严重问题，尝试修正
⋮----
confidence = review_result.get('confidence', 0.0)
⋮----
# 如果置信度高（>0.8），说明确实有问题
⋮----
# 尝试使用更多吸附位点重新生成
⋮----
# 这里可以实现重新生成逻辑
# 由于我们现在在可视化阶段，pathway已经生成完毕
# 最好的方式是记录问题，在下一次运行时使用更多位点
⋮----
# 记录问题到文件
issue_log = os.path.join(output_dir, f"surface_{miller_str}", "pathway", "structure_issues.txt")
⋮----
# 保存结构路径
⋮----
# 生成反应路径 GIF（按正确的反应序列）
⋮----
pathway_gif = self.viz_manager.visualize_reaction_pathway(
⋮----
sequence=reaction_intermediates,  # 保证顺序正确
⋮----
fps=1,  # 慢速播放以显示细节
hold_frames=5  # 每帧持续5帧以便观察
⋮----
# 5. 生成能量图 [05_energy_diagram]
⋮----
pathway = workflow_result.best_surface.pathway_analysis
⋮----
energies = []
labels = []
is_ts_list = []
⋮----
# 从pathway结果中提取完整能量剖面（包括中间体和过渡态）
# 对于每个反应步骤，添加：起始态 -> 过渡态 -> 产物态
⋮----
# 第一个step：添加起始态
⋮----
# 如果有能垒数据，添加过渡态
⋮----
ts_label = f"TS_{i+1}"
# TS能量 = 反应物能量 + 活化能
ts_energy = step.reactant_energy + step.activation_energy
⋮----
# 添加产物态（下一个step的起始态）
⋮----
# 提取barriers用于显示
barriers = {}
⋮----
energy_diagram = self.viz_manager.visualize_energy_diagram(
⋮----
# 6. 可视化KMC动力学轨迹（如果运行了KMC）[07_kmc_dynamics]
⋮----
# 加载中间体结构
⋮----
structures = {}
⋮----
# 从KMC结果创建轨迹
kmc_result = workflow_result.kmc_result
⋮----
state_history = kmc_result.state_history
time_history = kmc_result.time_history
⋮----
# 采样（最多100帧）
max_frames = 100
total_steps = len(state_history)
⋮----
indices = np.linspace(0, total_steps - 1, max_frames, dtype=int)
state_history_sampled = [state_history[i] for i in indices]
time_history_sampled = [time_history[i] for i in indices]
⋮----
state_history_sampled = state_history
time_history_sampled = time_history
⋮----
# 构建轨迹
trajectory = []
titles = []
⋮----
title = f"Frame {i+1}/{len(state_history_sampled)} | State: {state} | t = {time:.4e} s"
⋮----
# 渲染GIF
kmc_gif = os.path.join(output_dir, "visualizations", "07_kmc_surface_dynamics.gif")
⋮----
# 7. 生成HTML摘要 [06_html_summary] - 最后生成，包含所有可视化
⋮----
html_summary = self.viz_manager.generate_summary_html(
⋮----
# 打印摘要
⋮----
# =============================================================================
# 便捷函数
⋮----
"""
    便捷函数：运行带可视化的完整催化工作流

    Parameters
    ----------
    bulk_structure : str
        Bulk structure file path
    reaction_intermediates : List[str]
        Reaction intermediate list (e.g., ['*CO', '*O', '*CO2', '*'])
    initial_reactant : str
        Initial reactant (e.g., '*CO')
    temperature : float
        Temperature (K)
    **kwargs
        Additional parameters

    Returns
    -------
    Dict
        Results and visualization paths
    """
# 检查本地模型
fairchem_model = kwargs.get('fairchem_model', 'uma-s-1p1')
fairchem_model_path = kwargs.get('fairchem_model_path', None)
⋮----
# 如果没有提供本地路径，尝试使用项目中的模型
⋮----
project_model_path = os.path.join(
⋮----
fairchem_model_path = project_model_path
⋮----
dt = VisualizedGasSolidDigitalTwin(
⋮----
# 允许用户指定renderer，如果不指定则使用默认（tachyon）
````

## File: camel_agents/workflow.py
````python
"""CatDT CAMEL workflow: orchestration + runtime, helper logic delegated to tools."""
⋮----
_current_file_dir = Path(__file__).parent.resolve()
_project_root = _current_file_dir.parent
⋮----
logger = logging.getLogger(__name__)
⋮----
class CatDTCamelWorkflow
⋮----
"""CatDT CAMEL workflow orchestrator composed from dedicated modules."""
⋮----
def _make_memory(self, token_limit: int = 12000, window_size: int = 50) -> ChatHistoryMemory
⋮----
token_counter = getattr(self.llm, "token_counter", None)
⋮----
context_creator = ScoreBasedContextCreator(
⋮----
def _init_llm(self) -> None
⋮----
llm_temperature = float(os.getenv("OPENAI_TEMPERATURE", "0.2"))
llm_max_tokens = int(os.getenv("OPENAI_MAX_TOKENS", "4096"))
llm_top_p = float(os.getenv("OPENAI_TOP_P", "0.95"))
llm_presence_penalty = os.getenv("OPENAI_PRESENCE_PENALTY")
llm_frequency_penalty = os.getenv("OPENAI_FREQUENCY_PENALTY")
llm_seed = os.getenv("OPENAI_SEED")
⋮----
def _init_camel_tools(self) -> None
⋮----
registry = self.tools.to_camel_tools()
⋮----
required = {
missing = sorted(k for k in required if k not in registry)
⋮----
def _init_agents(self) -> None
⋮----
default_step_timeout = float(os.getenv("CATDT_AGENT_STEP_TIMEOUT_SEC", "600"))
agent45_step_timeout = float(os.getenv("CATDT_AGENT45_STEP_TIMEOUT_SEC", "1800"))
default_message_window = int(os.getenv("CATDT_AGENT_MESSAGE_WINDOW", "20"))
agent45_message_window = int(os.getenv("CATDT_AGENT45_MESSAGE_WINDOW", "8"))
default_token_limit = int(os.getenv("CATDT_AGENT_TOKEN_LIMIT", "200000"))
agent45_token_limit = int(os.getenv("CATDT_AGENT45_TOKEN_LIMIT", "200000"))
agent45_memory_token_limit = int(os.getenv("CATDT_AGENT45_MEMORY_TOKEN_LIMIT", "6000"))
agent45_memory_window = int(os.getenv("CATDT_AGENT45_MEMORY_WINDOW", "12"))
agent45_summarize_threshold = int(
agent45_memory_agent4 = self._make_memory(
agent45_memory_agent5 = self._make_memory(
⋮----
agent4_tools: List[FunctionTool] = [
agent5_tools: List[FunctionTool] = []  # Agent5 uses pre-check results only, no tool calls
⋮----
# Mechanism search agent (between Agent3 and Agent4/5)
mechanism_tools: List[FunctionTool] = [
⋮----
def _run_task(self, task: Task, use_memory: bool = False) -> TaskOutput
⋮----
_ = use_memory  # memory is managed by ChatAgent instances directly
allow_tool_calls = True if task.allow_tool_calls is None else bool(task.allow_tool_calls)
output: TaskOutput
retry_done = False
⋮----
output = task.agent.run(
⋮----
err_text = str(exc)
recoverable_tool_transcript_error = (
⋮----
retry_done = True
⋮----
chat_agent = getattr(task.agent, "chat_agent", None)
⋮----
pathway_dir = Path(state.output_base_dir) / state.run_id / "04_pathway"
prompts_dir = pathway_dir / "prompts"
⋮----
iter_file = prompts_dir / f"iter_{iteration:02d}_{agent_name}_prompt.txt"
⋮----
merged_file = pathway_dir / f"{agent_name}_prompts_full.txt"
⋮----
def _export_pathway_iteration_visuals(self, state: WorkflowState, iteration: int) -> None
⋮----
viz_dir = Path(state.output_base_dir) / state.run_id / "04_pathway" / f"iter_{iteration:02d}_step_visualizations"
⋮----
visualizer = CatalystSurfaceVisualizer(quality="low", auto_expand=False)
⋮----
step_name = str(step.get("name", f"step_{step_idx}"))
safe_name = re.sub(r"[^\w\-]+", "_", step_name)
reactant = step.get("reactant")
product = step.get("product")
⋮----
def _summarize_iteration_history(self, state: WorkflowState, max_items: int = 5) -> str
⋮----
lines: List[str] = []
⋮----
idx = item.get("iteration", "?")
status = item.get("status", "UNKNOWN")
issue_count = item.get("issue_count", 0)
issue_preview = item.get("issue_preview", "")
⋮----
@staticmethod
    def _truncate_prompt_block(text: Any, max_chars: int, block_name: str) -> str
⋮----
content = str(text or "").strip()
⋮----
keep = max(256, int(max_chars) - 128)
dropped = len(content) - keep
⋮----
def _reset_agent45_dialogue(self) -> None
⋮----
chat_agent = getattr(agent, "chat_agent", None)
⋮----
# === Mechanism search stage helpers ===
⋮----
@staticmethod
    def _should_enable_mechanism_stage(config) -> bool
⋮----
"""Check if the mechanism search stage should run."""
⋮----
"""Run mechanism search via Agent M1 (a real CAMEL agent with tools).

        The agent analyses the user's reaction description, decides the
        exploration mode, calls mechanism tools, and returns evaluated pathways.
        """
⋮----
# Build context for the agent prompt
initial_state = ""
target_state = ""
known_intermediates: List[str] = []
⋮----
intermediates = list(state.reaction_context.intermediates or [])
⋮----
initial_state = intermediates[0]
target_state = intermediates[-1]
known_intermediates = intermediates[1:-1]
⋮----
# Determine fixed site from current adsorbate position
fixed_site_str = "null"
⋮----
ads_atoms = read(state.reconstructed_surface_path)
ads_pos = ads_atoms.get_positions()[state.adsorbate_indices]
fixed_site = ads_pos[np.argmin(ads_pos[:, 2])].tolist()
fixed_site_str = str([round(x, 2) for x in fixed_site])
⋮----
surface_facet = ""
⋮----
m = _re.search(r"(\d{3,4})", str(state.surface_path))
⋮----
surface_facet = m.group(1)
⋮----
exploration_mode = getattr(config, "mechanism_exploration_mode", "agent_guided")
known_str = ", ".join(known_intermediates) if known_intermediates else "(none — recommend pathways)"
⋮----
prompt = f"""Analyse this catalytic reaction and run mechanism search.
⋮----
task = Task(
output = self._run_task(task)
⋮----
# Extract results from tool call records
# The tools have already stored results in self.tools._mechanism_context etc.
search_result = getattr(self.tools, "_mechanism_last_search_result", None)
shortlist = getattr(self.tools, "_mechanism_last_shortlist", None)
⋮----
# If agent didn't call all tools, fall back to direct tool calls
⋮----
"""Direct (non-agent) mechanism search pipeline as fallback."""
⋮----
runner = MechanismStageRunner(tools=self.tools, use_agents=False)
⋮----
"""Store mechanism search results into WorkflowState.

        Critically: updates ``state.reaction_context.intermediates`` with the
        best pathway so that Agent4/5 can consume it directly.
        """
⋮----
search_result = mechanism_output.get("search_result", {})
shortlist = mechanism_output.get("shortlist", {})
⋮----
all_pathways = search_result.get("pathways", [])
⋮----
# Update reaction_context.intermediates with the best pathway
# so Agent4/5 can consume it directly via existing _run_pathway_iteration
⋮----
best_pathway = state.retained_pathways[0]
best_intermediates = [
⋮----
def _build_agent45_transition_signature(self, state: WorkflowState) -> List[str]
⋮----
intermediates = self._get_working_intermediate_sequence(state)
⋮----
def _build_agent45_memento_query(self, state: WorkflowState, iteration: int) -> str
⋮----
reaction_type = str(state.reaction_context.reaction_type if state.reaction_context else "").strip()
constraints = str(state.reaction_context.constraints if state.reaction_context else "").strip()
transitions = self._build_agent45_transition_signature(state)
feedback = str(state.last_feedback or "").strip()
⋮----
feedback = feedback[:240] + "..."
query = (
# Add surface context so parametric retriever can match similar surfaces
⋮----
first_atoms = (state.step_structures[0] or {}).get("reactant")
⋮----
all_surf = list(state.surface_indices)
positions = first_atoms.get_positions()
syms = first_atoms.get_chemical_symbols()
surf_z = positions[all_surf, 2]
top_z = float(_np.max(surf_z))
top2_mask = surf_z >= (top_z - 4.0)
top2_indices = [all_surf[i] for i in range(len(all_surf)) if top2_mask[i]]
elem_counts = dict(Counter(syms[i] for i in top2_indices))
cell = first_atoms.cell
⋮----
def _retrieve_agent45_memento_context(self, state: WorkflowState, iteration: int) -> str
⋮----
query_text = self._build_agent45_memento_query(state=state, iteration=iteration)
⋮----
top_k = max(1, int(os.getenv("CATDT_AGENT45_MEMENTO_TOP_K", "6")))
min_score = float(os.getenv("CATDT_AGENT45_MEMENTO_MIN_SCORE", "0.05"))
include_negative = str(os.getenv("CATDT_AGENT45_MEMENTO_INCLUDE_NEGATIVE", "1")).strip().lower() not in {
⋮----
retrieval = self.tools.retrieve_agent45_memento_cases(
⋮----
def _retrieve_agent45_knowledge_context(self, state: WorkflowState, iteration: int) -> str
⋮----
query_text = state.memento_query_text or self._build_agent45_memento_query(state=state, iteration=iteration)
top_k = max(1, int(os.getenv("CATDT_AGENT45_KNOWLEDGE_TOP_K", "4")))
min_score = float(os.getenv("CATDT_AGENT45_KNOWLEDGE_MIN_SCORE", "0.03"))
⋮----
retrieval = self.tools.retrieve_agent45_knowledge_items(
⋮----
def _retrieve_agent45_skill_context(self, state: WorkflowState, iteration: int) -> str
⋮----
top_k = max(1, int(os.getenv("CATDT_AGENT45_SKILL_TOP_K", "4")))
min_score = float(os.getenv("CATDT_AGENT45_SKILL_MIN_SCORE", "0.03"))
⋮----
retrieval = self.tools.retrieve_agent45_skill_items(
⋮----
design_outline: List[Dict[str, Any]] = []
⋮----
step_outline: Dict[str, Any] = {
# Record actual positions Agent4 chose for added atoms
⋮----
# Build surface top-2-layer snapshot (same for all steps, stored once)
⋮----
surface_snapshot: Optional[Dict[str, Any]] = None
⋮----
# Use the first step's reactant to get surface coordinates
⋮----
surface_top_z = float(_np.max(surf_z))
# Top 2 layers: atoms within 4.0 Å of surface_top_z
top2_mask = surf_z >= (surface_top_z - 4.0)
⋮----
surface_snapshot = {
⋮----
# Build per-step adsorbate geometry snapshot from step_structures
# Stores core (original) and staged (Agent4-added) atoms separately
step_geometry_snapshot: List[Dict[str, Any]] = []
⋮----
snap: Dict[str, Any] = {"name": step.get("name", "")}
⋮----
atoms_obj = step.get(endpoint)
⋮----
ads_indices = list(step.get(f"{endpoint}_adsorbate_indices", []))
staged_indices = list(step.get(f"{endpoint}_staged_indices", []))
syms = atoms_obj.get_chemical_symbols()
pos = atoms_obj.get_positions()
staged_set = set(staged_indices)
core_ads = [
staged_atoms = []
core_positions = _np.array([a["xyz"] for a in core_ads]) if core_ads else _np.zeros((0, 3))
⋮----
sp = pos[si]
info: Dict[str, Any] = {
⋮----
dists = _np.linalg.norm(core_positions - sp, axis=1)
⋮----
neb_summary: Dict[str, Any] = {}
⋮----
def _strategy_hint(self) -> str
⋮----
hints = {
⋮----
def _create_reaction_context_task(self, state: WorkflowState) -> Task
⋮----
description = TaskPrompts.reaction_context_parsing(state.reaction_description)
⋮----
def _apply_reaction_context(output: TaskOutput)
⋮----
payload = output.json_dict or {}
context = output.pydantic
⋮----
context = ReactionContext(**payload)
⋮----
@staticmethod
    def _parse_pathway_step_formulas(step_name: str) -> Tuple[str, str]
⋮----
raw = str(step_name or "").strip()
⋮----
@staticmethod
    def _coerce_step_atoms_to_add(raw_step: Dict[str, Any]) -> List[Dict[str, Any]]
⋮----
normalized: List[Dict[str, Any]] = []
⋮----
raw_adds = raw_step.get("atoms_to_add", [])
⋮----
species = str(
⋮----
entry: Dict[str, Any] = {
⋮----
value = item.get(key)
⋮----
staged_atoms = raw_step.get("staged_atoms", [])
⋮----
position = item.get("position")
⋮----
pos = [float(position[0]), float(position[1]), float(position[2])]
⋮----
side = str(item.get("side") or "").strip().lower()
entry = {
⋮----
payload = dict(payload or {})
context_seq = (
⋮----
pathway_name = str(payload.get("pathway_name") or "").strip() or "repaired_pathway"
description_text = str(payload.get("description") or "").strip()
⋮----
summary = payload.get("summary")
summary_bits: List[str] = []
⋮----
validated = summary.get("validated_steps")
failed = summary.get("failed_steps")
⋮----
status = str(payload.get("status") or "").strip()
description_parts = [part for part in [pathway_name, status] if part]
⋮----
description_text = " | ".join(description_parts) or "Recovered Agent4 pathway design"
⋮----
overall_reaction = str(payload.get("overall_reaction") or "").strip()
⋮----
overall_reaction = " -> ".join(context_seq)
⋮----
overall_reaction = str(state.reaction_description or pathway_name)
⋮----
normalized_steps: List[Dict[str, Any]] = []
raw_steps = payload.get("steps", [])
⋮----
step_name = str(item.get("step_name") or item.get("name") or "").strip()
reactant_formula = str(item.get("reactant_formula") or "").strip()
product_formula = str(item.get("product_formula") or "").strip()
⋮----
reactant_formula = parsed_reactant
⋮----
product_formula = parsed_product
⋮----
reactant_formula = reactant_formula or str(context_seq[idx])
product_formula = product_formula or str(context_seq[idx + 1])
⋮----
step_name = f"{reactant_formula}_to_{product_formula}"
⋮----
confidence = item.get("confidence", 0.0)
⋮----
confidence_value = float(confidence or 0.0)
⋮----
confidence_value = 0.0
⋮----
validation = item.get("validation")
⋮----
interp = str(validation.get("interpolation") or "").strip().lower()
ok = bool(validation.get("ok"))
⋮----
confidence_value = 0.2
⋮----
confidence_value = 0.6
⋮----
confidence_value = 0.8
⋮----
raw_modify = item.get("atoms_to_modify", [])
atoms_to_modify = [
raw_remove = item.get("atoms_to_remove", [])
⋮----
confidence_value = payload.get("confidence", 0.0)
⋮----
top_confidence = float(confidence_value or 0.0)
⋮----
top_confidence = 0.0
⋮----
top_confidence = float(
⋮----
reaction_context = {
⋮----
# Surface info: include compact descriptors + full POSCAR text (LLM cannot read local files).
cell = structure.cell
surface_top_z = self._estimate_surface_top_z(structure, list(state.adsorbate_indices or []))
surface_element = structure[0].symbol if len(structure) > 0 else "unknown"
surface_poscar = self._atoms_to_poscar_text(structure)
surface_info = (
⋮----
# Add top-layer surface atom positions so LLM can find adsorption sites
⋮----
positions = structure.get_positions()
⋮----
# Top layer: atoms within 2.0 Å of surface_top_z
# (needs 2.0 Å because VSSR-MC reconstruction can add atoms
#  ~1.5 Å above the original slab top layer)
top_layer_mask = surf_z >= (surface_top_z - 2.0)
top_layer_indices = [all_surf[i] for i in range(len(all_surf)) if top_layer_mask[i]]
⋮----
# Per-step adsorbate coordinates (compact format from build_steps_payload)
baseline_steps = "(none)"
⋮----
baseline_steps = self.tools.build_steps_payload(self, state.tool_baseline_steps)
⋮----
# Build feedback: Agent5 text + previous round's adsorbate coords
feedback = state.last_feedback or ""
⋮----
feedback = (
⋮----
prev_adsorbate_summary = self.tools.build_steps_payload(self, state.step_structures)
feedback = f"{feedback}\n\n【前一轮各步吸附分子坐标】\n{prev_adsorbate_summary}"
⋮----
feedback_limit = int(os.getenv("CATDT_AGENT45_FEEDBACK_MAX_CHARS", "6000"))
history_limit = int(os.getenv("CATDT_AGENT45_HISTORY_MAX_CHARS", "1200"))
memento_limit = int(os.getenv("CATDT_AGENT45_MEMENTO_MAX_CHARS", "2000"))
knowledge_limit = int(os.getenv("CATDT_AGENT45_KNOWLEDGE_MAX_CHARS", "2000"))
skill_limit = int(os.getenv("CATDT_AGENT45_SKILL_MAX_CHARS", "2000"))
⋮----
feedback = self._truncate_prompt_block(feedback, feedback_limit, "agent4_feedback")
history_text = self._truncate_prompt_block(
memento_context = self._truncate_prompt_block(
knowledge_context = self._truncate_prompt_block(
skill_context = self._truncate_prompt_block(
⋮----
description = TaskPrompts.pathway_design(
⋮----
def _apply_pathway_design(output: TaskOutput)
⋮----
design = output.pydantic
⋮----
raw_text = str(output.raw or "").strip()
candidate = raw_text
fenced = re.search(r"```json\s*(\{.*?\})\s*```", raw_text, flags=re.DOTALL | re.IGNORECASE)
⋮----
candidate = fenced.group(1)
⋮----
loose = re.search(r"\{.*\}", raw_text, flags=re.DOTALL)
⋮----
candidate = loose.group(0)
⋮----
payload = json.loads(candidate)
⋮----
payload = payload or {}
⋮----
design = PathwayDesign(**self._coerce_pathway_design_payload(state, payload))
⋮----
repair_prompt = (
⋮----
repair_output = self.agent4.run(
repaired_design = repair_output.pydantic
⋮----
repair_payload = repair_output.json_dict or {}
⋮----
repair_raw = str(repair_output.raw or "").strip()
repair_candidate = repair_raw
repair_fenced = re.search(
⋮----
repair_candidate = repair_fenced.group(1)
⋮----
repair_loose = re.search(r"\{.*\}", repair_raw, flags=re.DOTALL)
⋮----
repair_candidate = repair_loose.group(0)
⋮----
repair_payload = json.loads(repair_candidate)
⋮----
repair_payload = repair_payload or {}
repaired_design = PathwayDesign(**self._coerce_pathway_design_payload(state, repair_payload))
design = repaired_design
⋮----
design = self._align_pathway_design_to_context(state, design)
⋮----
def _create_validation_task(self, state: WorkflowState) -> Task
⋮----
programmatic_check = self._agent5_programmatic_validation(state.step_structures)
energy_gate = dict(state.energy_gate_report or {})
design_steps = state.pathway_design.steps if state.pathway_design else []
structures_payload = self._build_agent5_structures_payload(state)
fatal_issues = list(programmatic_check.get("fatal_issues", []))
warning_issues = list(programmatic_check.get("warning_issues", []))
energy_fatal = list(energy_gate.get("fatal_issues", []))
energy_warning = list(energy_gate.get("warning_issues", []))
⋮----
precheck_summary = str(programmatic_check.get("summary", "(none)") or "(none)")
energy_summary = str(energy_gate.get("summary", "") or "").strip()
⋮----
precheck_summary = f"{precheck_summary}\n[UMA Energy Gate]\n{energy_summary}"
⋮----
step_briefs: List[str] = []
⋮----
spec = design_steps[idx - 1] if (idx - 1) < len(design_steps) else None
adds = len(spec.atoms_to_add) if spec else 0
rems = len(spec.atoms_to_remove) if spec else 0
⋮----
description = TaskPrompts.pathway_validation(
⋮----
def _apply_validation(output: TaskOutput)
⋮----
# LLM may return issues as List[dict] instead of List[str]; normalize
⋮----
normalized = []
⋮----
# Convert structured issue dict to readable string
parts = []
⋮----
report = output.pydantic
⋮----
report = ValidationReport(**payload)
⋮----
# Extra safety: normalize report.issues if pydantic auto-parsed but issues contain non-str
⋮----
combined_fatal_issues = list(fatal_issues)
combined_warning_issues = list(warning_issues)
⋮----
issues = list(report.issues or [])
⋮----
llm_status = str(report.status or "FAIL").upper()
llm_fail = llm_status != "PASS"
status = "FAIL" if combined_fatal_issues else "PASS"
⋮----
gate_note = str(programmatic_check.get("feedback", "")).strip() or (
⋮----
gate_note = f"{gate_note}\nUMA energy gate detected endpoint instability; revise staged coordinates."
feedback_text = str(report.feedback or "").strip()
⋮----
advisory = (
⋮----
issues = [f"[ADVISORY] {item}" for item in issues]
⋮----
# Prefix feedback with format reminder so Agent4 doesn't drift
⋮----
format_reminder = (
⋮----
def _validate_pathway_with_agent5(self, state: WorkflowState) -> Dict[str, Any]
⋮----
report = state.validation_report
⋮----
candidates: List[str] = []
⋮----
def _add(name: str) -> None
⋮----
reactant_raw = str(reactant or "").strip()
product_raw = str(product or "").strip()
reactant_clean = reactant_raw.lstrip("*")
product_clean = product_raw.lstrip("*")
⋮----
reactant_key = self._canonical_species_label(reactant_raw)
product_key = self._canonical_species_label(product_raw)
⋮----
struct_step_name = state.step_structures[step_index].get("name")
⋮----
neb_values = list(state.neb_results.values())
⋮----
@staticmethod
    def _extract_neb_activation_energy(neb_value: Any) -> Optional[float]
⋮----
value = None
⋮----
value = neb_value.get("activation_energy_forward", None)
⋮----
value = getattr(neb_value, "activation_energy_forward", None)
⋮----
f = float(value)
⋮----
def _collect_neb_quality_metrics(self, neb_results: Dict[str, Any]) -> Dict[str, float]
⋮----
total = len(neb_results or {})
⋮----
converged = 0
barriers: List[float] = []
⋮----
converged_flag = value.get("converged", None)
⋮----
converged_flag = getattr(value, "converged", None)
⋮----
ea = self._extract_neb_activation_energy(value)
⋮----
min_barrier = float(os.getenv("CATDT_AGENT45_BARRIER_MIN_EV", "0.05"))
max_barrier = float(os.getenv("CATDT_AGENT45_BARRIER_MAX_EV", "3.0"))
reasonable_count = sum(1 for ea in barriers if min_barrier <= ea <= max_barrier)
barrier_count = len(barriers)
⋮----
def _compute_iteration_reward(self, validation_passed: bool, neb_results: Dict[str, Any]) -> float
⋮----
metrics = self._collect_neb_quality_metrics(neb_results)
converged_ratio = float(metrics.get("converged_ratio", 0.0))
barrier_count = int(metrics.get("barrier_count", 0.0))
reasonable_ratio = float(metrics.get("reasonable_barrier_ratio", 0.0)) if barrier_count > 0 else 0.0
⋮----
reward = 0.4 * converged_ratio + 0.6 * reasonable_ratio
⋮----
def _run_pathway_iteration(self, state: WorkflowState, base_structure: Optional[Atoms], config: Optional[CatDTConfig]) -> None
⋮----
preplaced_products_cache: Optional[Dict[str, Dict[str, Any]]] = None
⋮----
iteration_num = i + 1
⋮----
preplaced_products_cache = self._agent4_tool_preplace_products(
preplaced_products = preplaced_products_cache
⋮----
# Expose baseline steps on tools so validate/suggest tools can access them
⋮----
baseline_dir = Path(state.output_base_dir) / state.run_id / "04_pathway" / f"iter_{iteration_num:02d}_tool_baseline"
⋮----
_ = self._retrieve_agent45_memento_context(state=state, iteration=iteration_num)
_ = self._retrieve_agent45_knowledge_context(state=state, iteration=iteration_num)
_ = self._retrieve_agent45_skill_context(state=state, iteration=iteration_num)
⋮----
agent4_retry_feedback = ""
max_agent4_attempts = 2
⋮----
design_task = self._create_pathway_design_task(
⋮----
design_issue = any(
⋮----
agent4_retry_feedback = (
⋮----
energy_gate_dir = Path(state.output_base_dir) / state.run_id / "04_pathway" / f"iter_{iteration_num:02d}_energy_gate"
⋮----
validation_task = self._create_validation_task(state)
⋮----
val_dir = Path(state.output_base_dir) / state.run_id / "04_pathway"
⋮----
val_payload = {
⋮----
validation_passed = bool(
⋮----
# Pre-NEB co-adsorbate relaxation: relax adsorbate+co-adsorbate with fixed surface
relax_dir = Path(state.output_base_dir) / state.run_id / "05_pre_neb_relax"
⋮----
neb_dir = Path(state.output_base_dir) / state.run_id / "06_neb"
⋮----
iter_reward = self._compute_iteration_reward(
⋮----
preview = ""
⋮----
preview = "; ".join(state.validation_report.issues[:2])
⋮----
"""Relax adsorbate atoms (surface fixed) on NEB endpoints before NEB.

        This ensures both endpoints are true local energy minima, which is critical
        for obtaining meaningful NEB barriers — especially when co-adsorbates have
        been placed at surface sites by Agent4.
        """
⋮----
predictor = self.tools.get_shared_fairchem_predictor(
⋮----
calc = predictor._calculator
⋮----
out_path = Path(output_dir)
⋮----
name = step.get("name", f"step{step_idx:02d}")
⋮----
atoms = step.get(side)
⋮----
ads_indices = set(step.get(f"{side}_adsorbate_indices", []))
staged_indices = set(step.get(f"{side}_staged_indices", []))
n_surface = int(step.get(f"{side}_fixed_atom_count", 0))
⋮----
# Identify top-2 surface layers by z-coordinate
top2_surface = set()
⋮----
surf_indices = [i for i in range(min(n_surface, len(atoms)))]
surf_z = np.array([atoms.positions[i][2] for i in surf_indices])
unique_z = np.sort(np.unique(np.round(surf_z, decimals=2)))
⋮----
z_threshold = unique_z[-2] - 0.1  # top 2 layers
⋮----
z_threshold = unique_z[-1] - 0.1  # only 1 layer
top2_surface = {i for i in surf_indices if atoms.positions[i][2] >= z_threshold}
⋮----
# Free atoms: adsorbate + staged + top-2 surface layers
free_indices = ads_indices | staged_indices | top2_surface
⋮----
original = atoms.copy()
# Track connectivity only for adsorbate + staged (not surface)
relax_free_for_connectivity = ads_indices | staged_indices
connectivity_before = self._free_atom_connectivity_signature(original, relax_free_for_connectivity)
a = atoms.copy()
⋮----
fixed_indices = [i for i in range(len(a)) if i not in free_indices]
⋮----
log_file = str(out_path / f"{name}_{side}_relax.log")
⋮----
e_before = a.get_potential_energy()
# Use a conservative optimizer step so staged atoms do not
# immediately collapse into new bonds before we can capture
# a lower-energy connectivity-preserving snapshot.
opt = BFGS(a, logfile=log_file)
last_safe = original.copy()
⋮----
safe_snapshot_steps = 0
connectivity_changed = False
⋮----
def _capture_last_safe_snapshot() -> None
⋮----
connectivity_now = self._free_atom_connectivity_signature(a, relax_free_for_connectivity)
⋮----
connectivity_changed = True
⋮----
last_safe = a.copy()
⋮----
safe_snapshot_steps = int(getattr(opt, "nsteps", 0))
⋮----
msg = str(exc)
⋮----
n_steps = int(getattr(opt, "nsteps", 0))
⋮----
selected = a
selected_connectivity = self._free_atom_connectivity_signature(selected, relax_free_for_connectivity)
⋮----
selected = last_safe.copy()
⋮----
chosen_label = (
⋮----
e_after = e_before
⋮----
selected_for_energy = selected.copy()
⋮----
e_after = selected_for_energy.get_potential_energy()
⋮----
# Update step structure in-place (remove constraint for NEB)
⋮----
valid = sorted({int(i) for i in free_indices if 0 <= int(i) < len(atoms)})
symbols = atoms.get_chemical_symbols()
signature: List[Tuple[int, int]] = []
⋮----
cutoff = 1.15 * (
⋮----
def _compute_episode_reward(self, state: WorkflowState) -> float
⋮----
validation_score = 1.0 if (state.validation_report and state.validation_report.status.upper() == "PASS") else 0.0
⋮----
neb_metrics = self._collect_neb_quality_metrics(state.neb_results or {})
neb_total = int(neb_metrics.get("total", 0.0))
⋮----
neb_score = (
⋮----
neb_score = 0.0
⋮----
target_steps = max(len((state.reaction_context.intermediates if state.reaction_context else [])) - 1, 1)
step_score = min(len(state.step_structures), target_steps) / target_steps
⋮----
reward = 0.5 * validation_score + 0.3 * neb_score + 0.2 * step_score
⋮----
# ==================================================================
# Multi-facet Wulff-weighted pipeline
⋮----
"""Run Agent3a→Agent2→Agent3b→Agent4/5→NEB→CatMAP for one facet.

        Mutates *state* in place. Returns a result dict with pathway_result_pkl,
        neb_summary, kmc_result_pkl, tof, and production_rates.
        """
is_multi = config.multi_facet
step_prefix = f"facet_{facet_id}/" if is_multi else ""
⋮----
# Initialize surface state for this facet
slab_atoms = read(surface_path)
⋮----
# ---- Agent2 output handler (closure over state) ----
def _apply_agent2_output(_output: TaskOutput)
⋮----
record = getattr(_output, "tool_call_record", None)
⋮----
tool_name = getattr(record, "tool_name", None)
args = getattr(record, "args", {}) or {}
result = getattr(record, "result", None)
⋮----
expected_surface = str(Path(state.surface_path).resolve())
actual_surface_raw = str(args.get("surface_path", "")).strip()
⋮----
actual_surface = str(Path(actual_surface_raw).resolve())
⋮----
ads_result_pkl = str(result or "").strip()
⋮----
ads_result = self.tools._load_result_from_pickle(ads_result_pkl)
best_adsorption_config_path = Path(ads_result_pkl).parent / "best_adsorption_config.vasp"
⋮----
ads_atoms = None
⋮----
best_result = getattr(ads_result, "best_result", None)
final_structure = getattr(best_result, "final_structure", None)
⋮----
ads_atoms = final_structure.copy()
⋮----
ads_atoms = read(state.surface_with_adsorbate_path)
⋮----
adsorbate_indices = list(range(min(state.surface_atom_count, len(ads_atoms)), len(ads_atoms)))
surface_indices = [i for i in range(len(ads_atoms)) if i not in set(adsorbate_indices)]
⋮----
# ---- Step 3a: Clean-slab reconstruction ----
⋮----
agent3a_sweeps = config.mc_clean_sweeps or max(1, config.mc_total_sweeps // 2)
⋮----
def _apply_agent3a_output(_output: TaskOutput)
⋮----
mc_result_pkl = str(result or "").strip()
⋮----
reconstructed_path = Path(mc_result_pkl).parent / "lowest_energy_reconstructed_structure.vasp"
recon_atoms = read(str(reconstructed_path))
⋮----
agent3a_task = Task(
⋮----
# ---- Step 2: Adsorption ----
agent2_task = Task(
⋮----
# ---- Step 3b: Reconstruction with adsorbate ----
def _apply_agent3b_output(_output: TaskOutput)
⋮----
expected_surface = str(Path(state.surface_with_adsorbate_path).resolve())
actual_surface_raw = str(args.get("surface_with_adsorbate_path", "")).strip()
⋮----
recon_atoms = read(state.reconstructed_surface_path)
reference_atoms = None
⋮----
reference_atoms = read(state.surface_with_adsorbate_path) if state.surface_with_adsorbate_path else None
⋮----
agent3b_step_name = "03b_reconstruction_adsorbate" if config.mc_two_step_reconstruction else "03_reconstruction"
agent3b_task = Task(
⋮----
base_structure = read(state.reconstructed_surface_path)
⋮----
# ---- Mechanism search (optional) ----
⋮----
mechanism_output = self._run_mechanism_search_stage(state, config)
⋮----
# ---- Agent4/5 pathway iteration ----
⋮----
# ---- NEB retry if needed ----
⋮----
neb_dir = Path(state.output_base_dir) / state.run_id / f"{step_prefix}06_neb"
⋮----
retry_results: Dict[str, Any] = {}
⋮----
step_result = self.tools.run_neb_for_steps(
⋮----
# ---- Build NEB summary ----
⋮----
# ---- Adsorption energies + pathway result ----
working_intermediates = self._get_working_intermediate_sequence(state)
adsorption_dir = Path(state.output_base_dir) / state.run_id / f"{step_prefix}04_adsorption_energies"
⋮----
fixed_site = None
⋮----
ads_pos = base_structure.get_positions()[state.adsorbate_indices]
⋮----
adsorption_intermediates: List[str] = []
⋮----
adsorption_intermediates = list(dict.fromkeys(working_intermediates))
⋮----
adsorption_result = self.tools.compute_adsorption_energies(
⋮----
transition_entries: List[Dict[str, str]] = []
⋮----
steps: List[PathwayStep] = []
⋮----
reactant = str(entry.get("reactant", "")).strip()
product = str(entry.get("product", "")).strip()
step_name = str(entry.get("step_name", f"{reactant}_to_{product}"))
reactant_energy = state.adsorbate_energies.get(reactant)
product_energy = state.adsorbate_energies.get(product)
reaction_energy = self._compute_reaction_energy_with_correction(
neb = self._find_neb_result_for_transition(state, reactant, product, i, step_name=step_name)
activation_energy = None
ts_energy = None
⋮----
activation_energy = getattr(neb, "activation_energy_forward", None)
⋮----
ts_energy = neb.energies[neb.transition_state_index]
⋮----
rate_determining_step = None
max_barrier = None
valid_barriers = [s for s in steps if s.activation_energy is not None]
⋮----
rate_determining_step = max(valid_barriers, key=lambda s: s.activation_energy)
⋮----
max_barrier = rate_determining_step.activation_energy
⋮----
overall_reaction_energy = 0.0
⋮----
first_e = state.adsorbate_energies.get(working_intermediates[0])
last_e = state.adsorbate_energies.get(working_intermediates[-1])
⋮----
overall_reaction_energy = last_e - first_e
⋮----
overall_reaction_energy = float(sum(s.reaction_energy for s in steps))
⋮----
pathway_dir = Path(state.output_base_dir) / state.run_id / f"{step_prefix}04_pathway"
⋮----
pathway_result = CompletePathwayResult(
pathway_pkl = self.tools._save_result_to_pickle(pathway_result, pathway_dir / "complete_pathway_result.pkl")
⋮----
# ---- KMC ----
kmc_result_pkl = None
tof = 0.0
production_rates: Dict[str, float] = {}
⋮----
kmc_result_pkl = self.tools.run_kmc_simulation(
kmc_result = self.tools._load_result_from_pickle(kmc_result_pkl)
tof = float(getattr(kmc_result, "tof", 0.0) or 0.0)
production_rates = dict(getattr(kmc_result, "production_rates", {}) or {})
⋮----
"""Run the per-facet pipeline for each top-N Wulff surface and aggregate."""
⋮----
miller_index = state.all_miller_indices.get(facet_id)
area_fraction = state.all_area_fractions.get(facet_id, 0.0)
⋮----
# Deep copy state so per-facet mutations are isolated
facet_state = copy.deepcopy(state)
# Clear per-facet mutable fields to avoid stale data
⋮----
result = self._run_single_facet_pipeline(
⋮----
@staticmethod
    def _aggregate_wulff_tof(state: WorkflowState) -> None
⋮----
"""Compute Wulff-weighted TOF: TOF_total = Σ (a_i / Σa_ok) × TOF_i."""
successful = [f for f in state.facet_results if f.error is None and f.tof > 0]
⋮----
total_area = sum(f.area_fraction for f in successful)
⋮----
all_species: set = set()
⋮----
def run(self, user_input_config: dict) -> dict
⋮----
config = CatDTConfig.create_from_args(**user_input_config)
⋮----
desired_output_dir = Path(config.output_base_dir).resolve()
current_output_dir = None
⋮----
current_output_dir = Path(getattr(self.tools, "output_base_dir", "")).resolve()
⋮----
# Sync electrochemical / model params into existing tools instance
⋮----
# Force re-creation of DT instance to pick up new params
⋮----
policy_path = Path(config.output_base_dir) / "evolvability_policy.json"
⋮----
state = WorkflowState(
⋮----
def _apply_agent1_output(_output: TaskOutput)
⋮----
slab_atoms = read(config.initial_surface_path)
⋮----
expected_bulk = str(Path(config.bulk_structure_path).resolve())
actual_bulk_raw = str(args.get("bulk_structure_path", "")).strip()
⋮----
actual_bulk = str(Path(actual_bulk_raw).resolve())
⋮----
surff_result_pkl = str(result or "").strip()
⋮----
surff_result = self.tools._load_result_from_pickle(surff_result_pkl)
top_surface = surff_result.get_top_n(1)[0]
miller_str = "".join(map(str, top_surface.miller_index))
⋮----
slab_atoms = read(state.surface_path)
⋮----
# Store all top-N surfaces for multi-facet mode
n_facets = config.multi_facet_top_n if config.multi_facet else 1
⋮----
ms = "".join(map(str, surf.miller_index))
slab_path = str(Path(surff_result_pkl).parent / f"slab_{ms}.vasp")
⋮----
agent1_task = Task(
⋮----
reaction_task = self._create_reaction_context_task(state)
⋮----
# --- Determine MC parameters (shared by all facets) ---
agent3_adsorbate_elements_for_mc = config.adsorbate_elements_for_mc
⋮----
slab = read(state.clean_slab_path)
symbols = slab.get_chemical_symbols()
⋮----
main_symbol = max(set(symbols), key=symbols.count)
agent3_adsorbate_elements_for_mc = [main_symbol]
agent3_mc_temperature = config.mc_temperature_k if config.mc_temperature_k is not None else 500.0
⋮----
# ================================================================
# Branch: multi-facet vs single-facet pipeline
⋮----
# Use the best facet's pathway for the final report
best_facet = max(
pathway_pkl = best_facet.pathway_result_pkl if best_facet else None
kmc_result_pkl = best_facet.kmc_result_pkl if best_facet else None
⋮----
results = {
⋮----
# Single-facet mode (original behavior)
facet_result = self._run_single_facet_pipeline(
⋮----
pathway_pkl = facet_result["pathway_result_pkl"]
kmc_result_pkl = facet_result.get("kmc_result_pkl")
⋮----
# Final report + policy update (shared by both modes)
⋮----
final_report = self.tools.generate_final_report(
⋮----
reward = self._compute_episode_reward(state)
⋮----
def __getattr__(self, name: str) -> Any
⋮----
tools = self.__dict__.get("tools")
⋮----
def main(user_input_config: dict)
⋮----
workflow = CatDTCamelWorkflow()
⋮----
__all__ = [
````

## File: core/active_learning/__init__.py
````python
__all__ = ["CandidateScore", "rank_candidates"]
````

## File: core/active_learning/acquisition.py
````python
def _minmax(values: np.ndarray) -> np.ndarray
⋮----
vmin = float(np.min(values))
vmax = float(np.max(values))
⋮----
def _pairwise_distance(a: np.ndarray, b: np.ndarray) -> float
⋮----
force = np.array([c.uncertainty_force for c in candidates], dtype=float)
energy = np.array([c.uncertainty_energy for c in candidates], dtype=float)
potential = np.array([c.uncertainty_potential for c in candidates], dtype=float)
⋮----
score_unc = (
⋮----
ranked_indices = np.argsort(-score_unc)
shortlist_indices = ranked_indices[: min(preselect_k, len(ranked_indices))]
shortlist = [candidates[int(i)] for i in shortlist_indices]
⋮----
selected: list[CandidateScore] = [shortlist[0]]
⋮----
best_idx = -1
best_score = -1.0
⋮----
unc_component = float(score_unc[int(shortlist_indices[idx])])
dist = min(
⋮----
_pairwise_distance(cand.descriptor, chosen.descriptor)  # type: ignore[arg-type]
⋮----
combined = (
⋮----
best_score = combined
best_idx = idx
````

## File: core/active_learning/framework.py
````python
@dataclass
class AdapterConfig
⋮----
type: str
params: dict[str, Any]
⋮----
@dataclass
class ALConfig
⋮----
workspace_dir: Path
train_file: Path
valid_file: Path
pool_file: Path
adapter: AdapterConfig
committee_seeds: list[int]
select_top_k: int
preselect_k: int
score_weights: dict[str, float]
diversity_lambda: float
⋮----
class BaseAdapter
⋮----
class CommandAdapter(BaseAdapter)
⋮----
def _run_shell(self, command: str, cwd: Path) -> None
⋮----
proc = subprocess.run(command, shell=True, cwd=str(cwd), check=False)
⋮----
train_tpl = cfg.adapter.params["train_command_template"]
ckpt_glob = cfg.adapter.params.get("checkpoint_glob", "*.model")
run_root = cfg.workspace_dir / "runs" / f"iter_{iteration}"
⋮----
checkpoints: list[Path] = []
⋮----
run_dir = run_root / f"seed_{seed}"
⋮----
command = train_tpl.format(
⋮----
models = sorted(
⋮----
score_tpl = cfg.adapter.params["score_command_template"]
score_file = cfg.workspace_dir / "artifacts" / "scores.json"
⋮----
command = score_tpl.format(
⋮----
payload = json.loads(score_file.read_text(encoding="utf-8"))
scores: list[CandidateScore] = []
⋮----
descriptor = (
⋮----
class CPMACEAdapter(BaseAdapter)
⋮----
command_adapter = CommandAdapter()
⋮----
cp_mace_root = Path(cfg.adapter.params["cp_mace_root"])
⋮----
device = cfg.adapter.params.get("device", "cuda")
dtype = cfg.adapter.params.get("default_dtype", "float64")
calc = MACECalculator(
⋮----
structures = ase.io.read(str(pool_file), index=":")
⋮----
_ = atoms.get_potential_energy()
⋮----
forces_comm = calc.results.get("forces_comm", None)
energy_var = calc.results.get("energy_var", 0.0)
⋮----
force_var = np.var(forces_comm, axis=0)
force_std = np.sqrt(np.sum(force_var, axis=1))
uncertainty_force = float(np.max(force_std))
⋮----
uncertainty_force = 0.0
uncertainty_energy = float(np.sqrt(max(float(energy_var), 0.0)))
⋮----
descriptor = None
⋮----
desc = calc.get_descriptors(atoms.copy(), invariants_only=True)
⋮----
structure_id = str(atoms.info.get("al_id", f"pool_{idx:08d}"))
⋮----
def _load_config(config_path: Path) -> ALConfig
⋮----
raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
workspace = Path(raw["workspace_dir"]).resolve()
⋮----
def _adapter_from_config(cfg: ALConfig) -> BaseAdapter
⋮----
def _read_pool(pool_file: Path) -> list[Any]
⋮----
rec_dir = cfg.workspace_dir / "recommendations"
⋮----
pool = _read_pool(pool_file)
⋮----
selected_id_set = {s.structure_id for s in selected}
by_index: dict[int, CandidateScore] = {}
⋮----
idx = int(s.metadata.get("index", -1))
⋮----
output_atoms = []
⋮----
sid = str(atoms.info.get("al_id", f"pool_{idx:08d}"))
⋮----
score = by_index.get(idx)
⋮----
rec_xyz = rec_dir / f"iter_{iteration:03d}_recommend.xyz"
⋮----
report = [
⋮----
def _save_state(cfg: ALConfig, state: dict[str, Any]) -> None
⋮----
def _load_state(cfg: ALConfig) -> dict[str, Any]
⋮----
state_file = cfg.workspace_dir / "state.json"
⋮----
def _append_xyz(base_xyz: Path, extra_xyz: Path, out_xyz: Path) -> None
⋮----
base_atoms = ase.io.read(str(base_xyz), index=":")
extra_atoms = ase.io.read(str(extra_xyz), index=":")
merged = list(base_atoms) + list(extra_atoms)
⋮----
pool_atoms = _read_pool(pool_file)
fb_atoms = ase.io.read(str(feedback_xyz), index=":")
fb_ids = {str(a.info["al_id"]) for a in fb_atoms if "al_id" in a.info}
⋮----
kept = []
⋮----
def initial_train_and_recommend(config_path: str) -> dict[str, Any]
⋮----
cfg = _load_config(Path(config_path).resolve())
adapter = _adapter_from_config(cfg)
⋮----
iteration = 0
checkpoints = adapter.train_committee(
scores = adapter.score_pool(cfg, checkpoints, cfg.pool_file)
selected = rank_candidates(
rec_xyz = _write_recommendations(cfg, iteration, selected, cfg.pool_file)
⋮----
state = {
⋮----
def retrain_with_feedback(config_path: str, feedback_xyz: str) -> dict[str, Any]
⋮----
state = _load_state(cfg)
⋮----
prev_iter = int(state["iteration"])
next_iter = prev_iter + 1
⋮----
train_prev = Path(state["train_file"]).resolve()
pool_prev = Path(state["pool_file"]).resolve()
feedback_file = Path(feedback_xyz).resolve()
⋮----
train_next = cfg.workspace_dir / "data" / f"train_iter_{next_iter:03d}.xyz"
pool_next = cfg.workspace_dir / "data" / f"pool_iter_{next_iter:03d}.xyz"
⋮----
checkpoints = adapter.train_committee(cfg, train_next, cfg.valid_file, next_iter)
scores = adapter.score_pool(cfg, checkpoints, pool_next)
⋮----
rec_xyz = _write_recommendations(cfg, next_iter, selected, pool_next)
⋮----
new_state = {
````

## File: core/active_learning/types.py
````python
@dataclass
class CandidateScore
⋮----
structure_id: str
uncertainty_force: float
uncertainty_energy: float
uncertainty_potential: float = 0.0
descriptor: Optional[np.ndarray] = None
metadata: dict[str, Any] = field(default_factory=dict)
⋮----
@property
    def uncertainty_total(self) -> float
````

## File: core/kmc/__init__.py
````python

````

## File: core/kmc/agent6_experiment_calibration.py
````python
def _rebuild_adsorbate_maps(pathway_result: CompletePathwayResult) -> None
⋮----
sequence_labels = []
sequence_energies = []
⋮----
adsorbate_energies = {
star_energy = adsorbate_energies.get("*", 0.0)
adsorption_energies = {
⋮----
def _update_rds(pathway_result: CompletePathwayResult) -> None
⋮----
valid = [step for step in pathway_result.steps if step.activation_energy is not None]
⋮----
reaction_energy_overrides = reaction_energy_overrides or {}
activation_energy_overrides = activation_energy_overrides or {}
⋮----
rebuilt = deepcopy(base_pathway_result)
⋮----
current_product_energy = float(base_pathway_result.steps[-1].product_energy)
⋮----
reaction_energy = float(reaction_energy_overrides.get(step.name, step.reaction_energy))
activation_energy = activation_energy_overrides.get(step.name, step.activation_energy)
activation_energy = None if activation_energy is None else float(activation_energy)
⋮----
current_product_energy = float(step.reactant_energy)
⋮----
candidates = [
````

## File: core/kmc/agent6_workflow_input_builder.py
````python
logger = logging.getLogger(__name__)
⋮----
_GAS_LABEL_RE = re.compile(r"^(?P<name>.+?)(?:\(g\)|_g)$")
_ITER_RE = re.compile(r"iter_(\d+)_energy_gate")
⋮----
@dataclass
class WorkflowStepRecord
⋮----
step_index: int
step_name: str
reactant_label: str
product_label: str
raw_reactant_energy: float
raw_product_energy: float
⋮----
def _is_gas_label(label: str) -> bool
⋮----
def _strip_gas_suffix(label: str) -> str
⋮----
match = _GAS_LABEL_RE.match(str(label).strip())
⋮----
def _normalize_surface_state_label(label: str) -> str
⋮----
label = str(label).strip()
⋮----
def _normalize_step_labels(reactant_raw: str, product_raw: str) -> Tuple[str, str]
⋮----
reactant_raw = str(reactant_raw).strip()
product_raw = str(product_raw).strip()
⋮----
reactant_is_gas = _is_gas_label(reactant_raw)
product_is_gas = _is_gas_label(product_raw)
⋮----
def _iteration_sort_key(path: Path) -> Tuple[int, str]
⋮----
match = _ITER_RE.search(path.parent.name)
⋮----
def _find_latest_energy_gate(run_dir: Path) -> Path
⋮----
candidates = sorted(
⋮----
def _load_step_records(energy_gate_path: Path) -> List[WorkflowStepRecord]
⋮----
payload = json.loads(energy_gate_path.read_text(encoding="utf-8"))
steps = payload.get("steps", [])
⋮----
records: List[WorkflowStepRecord] = []
⋮----
step_name = str(item.get("step_name", "")).strip()
⋮----
endpoints = item.get("endpoints", {})
reactant_energy = float(endpoints.get("reactant", {}).get("energy"))
product_energy = float(endpoints.get("product", {}).get("energy"))
⋮----
def _infer_surface_formula(run_dir: Path) -> Tuple[str, str]
⋮----
candidates = [
⋮----
atoms = read(str(candidate))
⋮----
fallback_candidates = list(run_dir.glob("05_pre_neb_relax/*_reactant_relaxed.vasp"))
⋮----
atoms = read(str(fallback_candidates[0]))
⋮----
def _extract_transition_state_from_traj(traj_path: Path) -> Dict[str, Any]
⋮----
frames = read(str(traj_path), index=":")
available: List[Tuple[int, float]] = []
⋮----
calc = getattr(atoms, "calc", None)
⋮----
results = getattr(calc, "results", {}) or {}
energy = results.get("energy")
⋮----
interior = [(idx, energy) for idx, energy in available if 0 < idx < len(frames) - 1]
source = interior if interior else available
⋮----
def _load_neb_endpoint_structures(traj_path: Path) -> Tuple[Any, Any]
⋮----
def _find_neb_traj(neb_dir: Path, step_name: str) -> Path
⋮----
candidate = neb_dir / f"{step_name}_neb.traj"
⋮----
def _build_adsorbate_energy_map(sequence_labels: Iterable[str], sequence_energies: Iterable[float]) -> Dict[str, float]
⋮----
energies: Dict[str, float] = {}
⋮----
def _build_adsorption_energy_map(adsorbate_energies: Dict[str, float]) -> Dict[str, float]
⋮----
star_energy = adsorbate_energies.get("*")
⋮----
run_path = Path(run_dir).resolve()
energy_gate_path = _find_latest_energy_gate(run_path)
neb_dir = run_path / "06_neb"
records = _load_step_records(energy_gate_path)
⋮----
step_diagnostics: List[Dict[str, Any]] = []
raw_step_entries: List[Dict[str, Any]] = []
⋮----
traj_path = _find_neb_traj(neb_dir, record.step_name)
⋮----
raw_reactant_energy = float(endpoint_energy_getter(record.step_name, "reactant", reactant_atoms))
raw_product_energy = float(endpoint_energy_getter(record.step_name, "product", product_atoms))
energy_source = "neb_endpoint_single_point"
⋮----
raw_reactant_energy = record.raw_reactant_energy
raw_product_energy = record.raw_product_energy
energy_source = "energy_gate"
⋮----
ts_info = _extract_transition_state_from_traj(traj_path)
raw_ts_energy = float(ts_info["transition_state_energy_raw"])
raw_barrier = max(0.0, raw_ts_energy - record.raw_reactant_energy)
⋮----
pathway_steps: List[PathwayStep] = []
sequence_labels: List[str] = []
sequence_energies: List[float] = []
⋮----
relative_product_energies: List[float] = [0.0] * len(raw_step_entries)
relative_reactant_energies: List[float] = [0.0] * len(raw_step_entries)
⋮----
current_product_energy = 0.0
⋮----
entry = raw_step_entries[rev_index]
⋮----
current_product_energy = relative_reactant_energies[rev_index]
⋮----
step_reactant_energy = relative_reactant_energies[index]
step_product_energy = relative_product_energies[index]
transition_state_energy = step_reactant_energy + float(entry["activation_energy_forward"])
⋮----
continuity_delta = 0.0
step_reactant_energy = float(entry["raw_reactant_energy"])
⋮----
previous_product_energy = pathway_steps[index - 1].product_energy
continuity_delta = abs(float(entry["raw_reactant_energy"]) - previous_product_energy)
step_reactant_energy = previous_product_energy
⋮----
step_product_energy = float(entry["raw_product_energy"])
⋮----
adsorbate_energies = _build_adsorbate_energy_map(sequence_labels, sequence_energies)
adsorption_energies = _build_adsorption_energy_map(adsorbate_energies)
⋮----
rate_determining_step = max(pathway_steps, key=lambda step: step.activation_energy or float("-inf"), default=None)
max_barrier = rate_determining_step.activation_energy if rate_determining_step is not None else None
⋮----
result = CompletePathwayResult(
⋮----
diagnostics = {
⋮----
out_dir = Path(output_dir).resolve()
⋮----
pathway_path = out_dir / "complete_pathway_result.pkl"
⋮----
summary = {
summary_path = out_dir / "workflow_agent6_input_summary.json"
⋮----
saved = {
⋮----
kmc_path = out_dir / "catmap_result.pkl"
⋮----
new_barrier = float(barrier_overrides[step.name])
⋮----
step_diag_map = {
⋮----
continuity_reactant_energy = float(step_diag_map[step_name].get("continuity_reactant_energy", 0.0))
⋮----
valid_barriers = [step for step in pathway_result.steps if step.activation_energy is not None]
⋮----
work_dir = str(Path(output_dir).resolve()) if output_dir else tempfile.mkdtemp(prefix="catdt-agent6-catmap-")
⋮----
predictor_factory = CatMAPPredictor
⋮----
dt = GasSolidDigitalTwin.__new__(GasSolidDigitalTwin)
⋮----
mechanism_builder = dt._build_catmap_mechanism
⋮----
predictor = predictor_factory(work_dir=work_dir)
⋮----
gas_names = list(getattr(predictor, "gases", {}).keys())
````

## File: core/kmc/catmap_predictor.py
````python
#!/usr/bin/env python
⋮----
"""
CatMAP Predictor - A wrapper for microkinetic modeling using CatMAP
"""
⋮----
# =============================================================================
# Data Classes for Input Specification
⋮----
@dataclass
class GasSpecies
⋮----
"""Gas phase species definition."""
name: str                          # e.g., 'CO', 'O2', 'H2O'
pressure: float = 1.0              # Partial pressure (bar)
formation_energy: float = 0.0      # Formation energy (eV), relative to references
frequencies: List[float] = field(default_factory=list)  # Vibrational frequencies (cm^-1)
⋮----
@dataclass
class AdsorbateSpecies
⋮----
"""Adsorbed species definition."""
name: str                          # e.g., 'CO', 'O', 'OH'
site: str = '111'                  # Adsorption site
formation_energies: Dict[str, float] = field(default_factory=dict)  # {surface: energy}
frequencies: List[float] = field(default_factory=list)
⋮----
@dataclass
class TransitionState
⋮----
"""Transition state definition."""
name: str                          # e.g., 'O-CO', 'H-OH'
site: str = '111'
⋮----
scaling_mode: str = 'initial_state'  # 'initial_state', 'final_state', 'BEP', or list
⋮----
@dataclass
class ReactionStep
⋮----
"""Elementary reaction step."""
expression: str                    # e.g., '*_s + CO_g -> CO*'
# Supports CatMAP syntax:
# - Adsorption: '*_s + CO_g -> CO*'
# - Dissociation: '2*_s + O2_g <-> O-O* + *_s -> 2O*'
# - Surface reaction: 'CO* + O* <-> O-CO* + * -> CO2_g + 2*'
# - Electrochemical: 'H* + pe_g -> H2_g + *_s' (pe = proton-electron pair)
⋮----
# Main Predictor Class
⋮----
class CatMAPPredictor
⋮----
"""
    Wrapper class for CatMAP microkinetic modeling.

    Simplifies the process of setting up and running microkinetic models
    for arbitrary reaction systems.

    Example usage:
    -------------
    ```python
    predictor = CatMAPPredictor()

    # Define reaction mechanism
    predictor.add_reaction('*_s + CO_g -> CO*')
    predictor.add_reaction('2*_s + O2_g <-> O-O* + *_s -> 2O*')
    predictor.add_reaction('CO* + O* <-> O-CO* + * -> CO2_g + 2*')

    # Add energy data
    predictor.add_gas('CO', formation_energy=0.0, frequencies=[2170])
    predictor.add_gas('O2', formation_energy=0.0, frequencies=[1580])
    predictor.add_gas('CO2', formation_energy=-2.45, frequencies=[1333, 2349, 667, 667])

    predictor.add_adsorbate('CO', {'Pt': 1.7, 'Pd': 1.55, 'Cu': 2.58})
    predictor.add_adsorbate('O', {'Pt': 1.62, 'Pd': 1.55, 'Cu': 1.07})
    predictor.add_transition_state('O-CO', {'Pt': 4.04, 'Pd': 4.2, 'Cu': 4.18})

    # Set conditions and run
    predictor.set_conditions(temperature=500, pressures={'CO_g': 1.0, 'O2_g': 0.33})
    result = predictor.run()
    ```
    """
⋮----
# Available thermodynamic modes
GAS_THERMO_MODES = [
⋮----
'shomate_gas',      # Shomate polynomials (most accurate)
'ideal_gas',        # Ideal gas with statistical mechanics
'zero_point_gas',   # Only zero-point corrections
'fixed_entropy_gas', # Fixed entropy approximation
'frozen_gas',       # No thermal contributions
⋮----
ADSORBATE_THERMO_MODES = [
⋮----
'harmonic_adsorbate',   # Harmonic oscillator
'hindered_adsorbate',   # Hindered translator/rotor
'zero_point_adsorbate', # Only zero-point corrections
'frozen_adsorbate',     # No thermal contributions
⋮----
ELECTROCHEMICAL_THERMO_MODES = [
⋮----
'simple_electrochemical',    # Basic CHE model
'hbond_electrochemical',     # With H-bonding corrections
⋮----
'field_electrochemical',     # With field effects
⋮----
# Available output variables
OUTPUT_VARIABLES = [
⋮----
# Solver level outputs
'coverage',           # Surface coverages
'rate',               # Reaction rates
'production_rate',    # Production rates of species
'consumption_rate',   # Consumption rates
'turnover_frequency', # TOF
'selectivity',        # Product selectivity
'carbon_selectivity', # Carbon-based selectivity
'rate_control',       # Degree of rate control
'selectivity_control',# Degree of selectivity control
'rxn_direction',      # Net reaction direction
'directional_rates',  # Forward/reverse rates
'rxn_order',          # Apparent reaction orders
⋮----
# Scaler level outputs
'rxn_parameter',      # Reaction parameters
'frequency',          # Frequencies
'electronic_energy',  # Electronic energies
'free_energy',        # Free energies
'zero_point_energy',  # ZPE
'enthalpy',           # Enthalpies
'entropy',            # Entropies
⋮----
def __init__(self, catmap_root: Optional[str] = None, work_dir: Optional[str] = None)
⋮----
"""
        Initialize CatMAP predictor.

        Args:
            catmap_root: Path to CatMAP installation (auto-detected if None)
            work_dir: Working directory for output files
        """
⋮----
# Storage for model components
⋮----
# Model settings
self.temperature: float = 500.0  # K
self.voltage: Optional[float] = None  # V vs SHE (for electrochemistry)
self.beta: float = 0.5  # Symmetry factor
⋮----
# Descriptor settings
⋮----
# Thermodynamic settings
⋮----
# Scaling settings
⋮----
# Interaction settings
self.adsorbate_interaction_model: Optional[str] = None  # None, 'ideal', 'first_order'
⋮----
# Output settings
⋮----
# Site definition
⋮----
# Atomic reservoir settings (for defining references for each element)
⋮----
# Solver settings
⋮----
# Results storage
⋮----
# =========================================================================
# Methods for adding model components
⋮----
def add_reaction(self, expression: str)
⋮----
"""
        Add an elementary reaction step.

        Supports CatMAP syntax:
        - Adsorption: '*_s + CO_g -> CO*' or '*_s + CO_g -> CO_s'
        - Desorption: 'CO* -> *_s + CO_g'
        - Dissociation with TS: '2*_s + O2_g <-> O-O* + *_s -> 2O*'
        - Surface reaction: 'CO* + O* <-> O-CO* + * -> CO2_g + 2*'
        - Electrochemical: 'H* + pe_g -> H2_g + *_s' (pe = proton-electron)

        Args:
            expression: Reaction expression string
        """
⋮----
"""
        Add a gas phase species.

        Args:
            name: Species name (e.g., 'CO', 'O2', 'H2O')
            formation_energy: Formation energy in eV (relative to references)
            pressure: Partial pressure in bar
            frequencies: Vibrational frequencies in cm^-1
        """
⋮----
"""
        Add an adsorbed species.

        Args:
            name: Adsorbate name (e.g., 'CO', 'O', 'OH')
            formation_energies: Dict mapping surface names to formation energies (eV)
            site: Adsorption site type (default: '111')
            frequencies: Vibrational frequencies in cm^-1
        """
⋮----
# Auto-add surfaces
⋮----
"""
        Add a transition state.

        Args:
            name: TS name (e.g., 'O-CO', 'H-OH') - use '-' to indicate bond forming/breaking
            formation_energies: Dict mapping surface names to formation energies (eV)
            site: Site type
            frequencies: Vibrational frequencies
            scaling_mode: 'initial_state', 'final_state', 'BEP', or explicit scaling list
        """
⋮----
# Methods for setting conditions
⋮----
"""
        Set reaction conditions.

        Args:
            temperature: Temperature in K
            pressures: Dict of gas partial pressures {species_g: pressure}
            voltage: Electrode potential in V vs SHE (for electrochemistry)
            pH: Solution pH (for electrochemistry)
            beta: Symmetry factor for electrochemical reactions
        """
⋮----
# Remove '_g' suffix if present
name = species.replace('_g', '')
⋮----
def set_atomic_reservoirs(self, reservoir_list: List[str])
⋮----
"""
        Set the atomic reservoir list for CatMAP.

        This tells CatMAP which gas-phase species to use as references for
        each element in the reaction mechanism.

        Args:
            reservoir_list: List of gas species (e.g., ['H2_g', 'CO_g', 'H2O_g'])
        """
⋮----
"""
        Set descriptor space for volcano plots.

        Args:
            names: Descriptor names (must match adsorbate names, e.g., ['O_s', 'CO_s'])
            ranges: Descriptor ranges [[min1, max1], [min2, max2]]
            resolution: Grid resolution
        """
⋮----
"""
        Set thermodynamic correction modes.

        Args:
            gas_mode: One of GAS_THERMO_MODES
            adsorbate_mode: One of ADSORBATE_THERMO_MODES
            electrochemical_mode: One of ELECTROCHEMICAL_THERMO_MODES
        """
⋮----
def set_scaling_constraints(self, constraints: Dict[str, Union[str, List]])
⋮----
"""
        Set scaling constraints for adsorbates and transition states.

        Args:
            constraints: Dict mapping species to constraints
                - For adsorbates: [slope1, slope2, intercept] or ['+', 0, None]
                  '+' means fit slope, 0 means fixed at 0, None means fit intercept
                - For TS: 'initial_state', 'final_state', 'BEP', or explicit list

        Example:
            predictor.set_scaling_constraints({
                'O_s': ['+', 0, None],      # Scales only with descriptor 1
                'CO_s': [0, '+', None],     # Scales only with descriptor 2
                'O-CO_s': 'initial_state',  # TS scales with IS
            })
        """
⋮----
"""
        Set adsorbate-adsorbate interaction parameters.

        Args:
            model: 'first_order' or None
            self_interactions: {species: [param_per_surface]}
            cross_interactions: {species1: {species2: [param_per_surface]}}
        """
⋮----
# Input file generation
⋮----
def _generate_energies_file(self, filepath: str)
⋮----
"""Generate the energies.txt input file."""
lines = ['surface_name\tsite_name\tspecies_name\tformation_energy\tbulk_structure\tfrequencies\tother_parameters\treference']
⋮----
# Add gas species
⋮----
freq_str = str(gas.frequencies) if gas.frequencies else '[]'
⋮----
# Add adsorbates
⋮----
freq_str = str(ads.frequencies) if ads.frequencies else '[]'
⋮----
# Add transition states
⋮----
freq_str = str(ts.frequencies) if ts.frequencies else '[]'
⋮----
def _generate_mkm_file(self, filepath: str, energies_file: str)
⋮----
"""Generate the .mkm setup file."""
lines = []
⋮----
# Reaction expressions
⋮----
# Surface names
⋮----
# Descriptor names and ranges
⋮----
# Auto-detect descriptors from adsorbates (use first two)
ads_names = list(self.adsorbates.keys())[:2]
desc_names = [f"{name}_s" for name in ads_names]
⋮----
# Auto-determine ranges from formation energies
ranges = []
⋮----
energies = list(self.adsorbates[name].formation_energies.values())
⋮----
# Temperature
⋮----
# Electrochemical settings
⋮----
# Species definitions
⋮----
# Gas pressures
⋮----
# Add proton-electron if electrochemical
⋮----
# Files
⋮----
# Thermodynamic modes
⋮----
# Scaling constraints
⋮----
# Auto-generate default constraints based on descriptor order
⋮----
# Get descriptor base names (without _s)
desc_base_names = [d.replace('_s', '') for d in self.descriptor_names] if self.descriptor_names else []
⋮----
constraint = [0] * len(desc_base_names) + [None]
⋮----
idx = desc_base_names.index(name)
⋮----
mode = ts.scaling_mode
⋮----
# Interaction model
⋮----
# Atomic reservoir list (if set) - used in _generate_mkm_file
⋮----
# Running the model
⋮----
def run(self, output_dir: Optional[str] = None) -> 'CatMAPResult'
⋮----
"""
        Run the microkinetic model.

        Args:
            output_dir: Directory for output files (uses work_dir if None)

        Returns:
            CatMAPResult object containing all results
        """
⋮----
output_dir = output_dir or self.work_dir
⋮----
# Generate input files
energies_file = os.path.join(output_dir, 'energies.txt')
mkm_file = os.path.join(output_dir, 'model.mkm')
⋮----
# Change to output directory and run
original_dir = os.getcwd()
⋮----
# Create and run model
⋮----
# Create result object
⋮----
def run_single_point(self, surface: str, output_dir: Optional[str] = None) -> 'SinglePointResult'
⋮----
"""
        Run microkinetic model for a single surface (single point calculation).

        This is useful when you want to calculate rates/coverages for a specific
        catalyst without scanning the entire descriptor space.

        Args:
            surface: Surface name (must have energies defined)
            output_dir: Directory for output files

        Returns:
            SinglePointResult object with rates and coverages
        """
⋮----
mkm_file = os.path.join(output_dir, 'single_point.mkm')
⋮----
model = ReactionModel(setup_file='single_point.mkm')
⋮----
# Extract results at the single point
result = SinglePointResult(
⋮----
def _generate_single_point_mkm(self, filepath: str, energies_file: str, surface: str)
⋮----
"""Generate mkm file for single point calculation."""
⋮----
# Only the specified surface
⋮----
# Get descriptors for this surface from adsorbate energies
desc_values = []
⋮----
desc_base = desc.replace('_s', '')
⋮----
energy = self.adsorbates[desc_base].formation_energies.get(surface)
⋮----
# Auto-detect from first two adsorbates
ads_list = list(self.adsorbates.keys())[:2]
⋮----
desc_names = []
⋮----
energy = self.adsorbates[ads].formation_energies.get(surface)
⋮----
# Single point: resolution=1, range is just the descriptor value
⋮----
ranges = [[v, v] for v in desc_values]
⋮----
# Scaling - for single surface, use direct values
⋮----
desc_base_names = [d.replace('_s', '') for d in self.descriptor_names]
⋮----
# Atomic reservoir list (if set)
⋮----
# Single Point Result Class
⋮----
class SinglePointResult
⋮----
"""Container for single point calculation results."""
⋮----
def __init__(self, surface: str, model, output_dir: str, predictor)
⋮----
# Extract results
⋮----
def __getstate__(self)
⋮----
# ReactionModel and predictor runtime objects are not stable pickle payloads.
# Persist only the extracted single-point observables plus minimal metadata.
⋮----
def __setstate__(self, state)
⋮----
def _extract_results(self)
⋮----
"""Extract coverages and rates from the model."""
# Get coverage at the single point
coverage_map = getattr(self.model, 'coverage_map', None)
rate_map = getattr(self.model, 'rate_map', None)
production_rate_map = getattr(self.model, 'production_rate_map', None)
⋮----
# coverage_map is list of [descriptors, [coverages]]
point = coverage_map[0]
coverage_values = point[1] if len(point) > 1 else []
⋮----
# Map to species names
adsorbate_names = self.model.adsorbate_names
⋮----
point = rate_map[0]
rate_values = point[1] if len(point) > 1 else []
⋮----
# Map to reaction names
rxn_expressions = self.model.rxn_expressions
⋮----
point = production_rate_map[0]
prod_values = point[1] if len(point) > 1 else []
⋮----
# Map to gas species
gas_names = self.model.gas_names
⋮----
def get_coverage(self, species: Optional[str] = None) -> Union[Dict, float]
⋮----
"""
        Get surface coverage.

        Args:
            species: Species name (e.g., 'CO_s'). None returns all.

        Returns:
            Coverage value or dict of all coverages
        """
⋮----
def get_rate(self, reaction: Optional[int] = None) -> Union[Dict, float]
⋮----
"""
        Get reaction rate.

        Args:
            reaction: Reaction index (1-indexed). None returns all.

        Returns:
            Rate value or dict of all rates
        """
⋮----
def get_production_rate(self, species: Optional[str] = None) -> Union[Dict, float]
⋮----
"""
        Get production rate.

        Args:
            species: Gas species name (e.g., 'CO2_g'). None returns all.

        Returns:
            Production rate or dict of all rates
        """
⋮----
def get_turnover_frequency(self) -> float
⋮----
"""Get the overall turnover frequency (max production rate)."""
⋮----
def summary(self) -> str
⋮----
"""Generate a text summary of results."""
lines = ["=" * 60]
⋮----
# Result Class
⋮----
class CatMAPResult
⋮----
"""Container for CatMAP results with analysis methods."""
⋮----
def __init__(self, model, output_dir: str, predictor: CatMAPPredictor)
⋮----
def get_coverage_map(self)
⋮----
"""Get coverage map data."""
⋮----
def get_rate_map(self)
⋮----
"""Get rate map data."""
⋮----
def get_production_rate_map(self)
⋮----
"""Get production rate map data."""
⋮----
def get_descriptor_dict(self)
⋮----
"""Get descriptor dictionary mapping surfaces to descriptor values."""
⋮----
def plot_rate(self, save: str = 'rate.png', **kwargs)
⋮----
"""Plot reaction rate volcano."""
⋮----
vm = analyze.VectorMap(self.model)
⋮----
save_path = os.path.join(self.output_dir, save)
⋮----
def plot_coverage(self, save: str = 'coverage.png', species: Optional[List[str]] = None, **kwargs)
⋮----
"""
        Plot surface coverage.

        Args:
            save: Output filename
            species: List of species to include (e.g., ['CO_s', 'O_s']). None for all.
        """
⋮----
def plot_production_rate(self, save: str = 'production_rate.png', **kwargs)
⋮----
"""Plot production rate."""
⋮----
def plot_scaling(self, save: str = 'scaling.png')
⋮----
"""Plot scaling relations."""
⋮----
sa = analyze.ScalingAnalysis(self.model)
⋮----
"""
        Plot potential energy diagram for reaction mechanism.

        Args:
            surfaces: List of surfaces to plot (default: all)
            mechanism: List of reaction step indices (1-indexed, negative for reverse)
            save: Output filename
        """
⋮----
ma = analyze.MechanismAnalysis(self.model)
⋮----
# Configure
⋮----
def plot_rate_control(self, save: str = 'rate_control.png', **kwargs)
⋮----
"""Plot degree of rate control."""
⋮----
# Need to recalculate with rate_control output
⋮----
mm = analyze.MatrixMap(self.model)
⋮----
def plot_selectivity(self, save: str = 'selectivity.png', **kwargs)
⋮----
"""Plot product selectivity."""
⋮----
def plot_all(self, prefix: str = '')
⋮----
"""Generate all standard plots."""
plots = {}
⋮----
def get_rate_at_surface(self, surface: str) -> Dict
⋮----
"""
        Get reaction rates at a specific surface.

        Args:
            surface: Surface name (e.g., 'Pt')

        Returns:
            Dict with rates for each reaction
        """
desc_dict = self.get_descriptor_dict()
⋮----
descriptors = desc_dict[surface]
# Find nearest point in rate_map
rate_map = self.get_rate_map()
⋮----
# rate_map is list of [descriptors, rates]
min_dist = float('inf')
best_rates = None
⋮----
dist = sum((d1 - d2)**2 for d1, d2 in zip(desc, descriptors))
⋮----
min_dist = dist
best_rates = rates
⋮----
# Convenience functions for common reactions
⋮----
"""
    Create a CO oxidation microkinetic model.

    Args:
        energies: Dict with keys 'CO', 'O', 'O-CO', 'O-O' mapping to {surface: energy}
        temperature: Reaction temperature in K
        co_pressure: CO partial pressure
        o2_pressure: O2 partial pressure

    Returns:
        Configured CatMAPPredictor
    """
predictor = CatMAPPredictor()
⋮----
# Add reactions
⋮----
# Add gases
⋮----
# Add adsorbates and TS
⋮----
"""
    Create a HER (Hydrogen Evolution Reaction) microkinetic model.

    Args:
        energies: Dict with keys 'H', 'H-H', 'pe-H' mapping to {surface: energy}
        temperature: Reaction temperature in K
        voltage: Electrode potential vs SHE
        pH: Solution pH

    Returns:
        Configured CatMAPPredictor
    """
⋮----
# Add reactions (Volmer-Heyrovsky-Tafel mechanism)
predictor.add_reaction('pe_g + *_s <-> H*')                           # Volmer
predictor.add_reaction('pe_g + H* <-> pe-H* <-> H2_g + *_s')          # Heyrovsky
predictor.add_reaction('H* + H* <-> H-H* + *_s <-> H2_g + 2*_s')      # Tafel
⋮----
"""
    Create an ORR (Oxygen Reduction Reaction) microkinetic model.

    Supports two mechanisms:
    - 'associative': O2 -> OOH -> O + H2O -> OH -> H2O (4e- pathway)
    - 'dissociative': O2 -> 2O -> 2OH -> 2H2O

    Args:
        energies: Dict with keys for intermediates mapping to {surface: energy}
            Required keys depend on mechanism:
            - associative: 'OH', 'O', 'OOH' (and optionally TS)
            - dissociative: 'OH', 'O', 'O-O'
        temperature: Reaction temperature in K
        voltage: Electrode potential vs SHE
        pH: Solution pH
        mechanism: 'associative' or 'dissociative'

    Returns:
        Configured CatMAPPredictor
    """
⋮----
# 4-electron associative pathway (most common on Pt-group metals)
# O2 + * + H+ + e- -> OOH*
# OOH* + H+ + e- -> O* + H2O
# O* + H+ + e- -> OH*
# OH* + H+ + e- -> H2O + *
⋮----
# Set descriptors (typically OH and O binding)
⋮----
else:  # dissociative
# Dissociative pathway
# O2 + 2* -> 2O*
⋮----
"""
    Create a CO2RR (CO2 Reduction Reaction) microkinetic model.

    Supports different product pathways:
    - 'CO': CO2 -> COOH -> CO (2e- pathway)
    - 'HCOOH': CO2 -> OCHO -> HCOOH (2e- pathway)
    - 'CH4': Full 8e- pathway to methane (simplified)

    Args:
        energies: Dict with keys for intermediates mapping to {surface: energy}
            Required keys depend on product:
            - CO: 'COOH', 'CO'
            - HCOOH: 'OCHO', 'HCOOH' (or use 'OCHO' only)
            - CH4: 'COOH', 'CO', 'CHO', 'CH2O', 'CH3O', 'O', 'OH'
        temperature: Reaction temperature in K
        voltage: Electrode potential vs SHE
        pH: Solution pH
        product: 'CO', 'HCOOH', or 'CH4'

    Returns:
        Configured CatMAPPredictor
    """
⋮----
# 2e- pathway to CO
# CO2 + * + H+ + e- -> COOH*
# COOH* + H+ + e- -> CO* + H2O
# CO* -> CO + *
⋮----
# 2e- pathway to formic acid/formate
# CO2 + * + H+ + e- -> OCHO*
# OCHO* + H+ + e- -> HCOOH + *
⋮----
# For single adsorbate, use it as both descriptors or add H
⋮----
# Simplified 8e- pathway to methane (via CO)
# CO2 -> COOH -> CO -> CHO -> CH2O -> CH3O -> CH4 + O -> OH -> H2O
⋮----
"""
    Create an NRR (Nitrogen Reduction Reaction) microkinetic model.

    Supports two mechanisms:
    - 'associative': Distal or alternating hydrogenation
    - 'dissociative': N2 dissociation followed by hydrogenation

    Args:
        energies: Dict with keys for intermediates mapping to {surface: energy}
            Typical keys: 'N2H', 'N2H2', 'N', 'NH', 'NH2', 'NH3'
        temperature: Reaction temperature in K
        voltage: Electrode potential vs SHE
        pH: Solution pH
        mechanism: 'associative' or 'dissociative'

    Returns:
        Configured CatMAPPredictor
    """
⋮----
# N2 + 2* -> 2N*
# N* + H+ + e- -> NH*
# NH* + H+ + e- -> NH2*
# NH2* + H+ + e- -> NH3 + *
⋮----
else:  # associative (distal)
# Associative distal pathway
# N2 + * -> N2*
# N2* + H+ + e- -> N2H*
# N2H* + H+ + e- -> N + NH3 (or further intermediates)
⋮----
# Main entry point for testing
⋮----
# Example: CO oxidation on various metal surfaces
# Data from Falsig et al. and Angew. Chem. Int. Ed., 47, 4835 (2008)
energies = {
⋮----
predictor = create_co_oxidation_model(energies, temperature=500)
⋮----
# Set descriptors explicitly
⋮----
# Test 1: Full descriptor space scan (volcano plot)
⋮----
result = predictor.run(output_dir='output/catmap_test')
⋮----
# Generate plots
⋮----
# Test 2: Single point calculation for Pt
⋮----
single_result = predictor.run_single_point('Pt', output_dir='output/catmap_single_pt')
⋮----
# Test 3: Compare multiple surfaces
⋮----
sp_result = predictor.run_single_point(surface, output_dir=f'output/catmap_single_{surface.lower()}')
tof = sp_result.get_turnover_frequency()
````

## File: core/kmc/pmutt_predictor.py
````python
#!/usr/bin/env python
⋮----
"""
pMuTT Predictor - A simplified wrapper for multiscale thermochemistry calculations
"""
⋮----
# =============================================================================
# Data Classes for Input Specification
⋮----
@dataclass
class SpeciesData
⋮----
"""Species thermodynamic data from DFT calculations."""
name: str                                  # Species name
atoms: Optional[object] = None             # ASE Atoms object
potentialenergy: Optional[float] = None    # DFT energy (eV)
vib_wavenumbers: List[float] = field(default_factory=list)  # cm^-1
symmetrynumber: int = 1                    # Symmetry number
spin: float = 0                            # Electronic spin
phase: str = 'gas'                         # 'gas' or 'surface'
elements: Dict[str, int] = field(default_factory=dict)  # {'H': 2, 'O': 1}
geometry: str = 'nonlinear'                # 'nonlinear', 'linear', 'atom'
# Optional coverage effects
coverage_effects: List[Dict] = field(default_factory=list)
# Optional notes
notes: str = ''
⋮----
@dataclass
class ReactionData
⋮----
"""Reaction data."""
reactants: List[str]                       # Reactant names
products: List[str]                        # Product names
reactants_stoich: List[float] = field(default_factory=list)
products_stoich: List[float] = field(default_factory=list)
transition_state: Optional[str] = None     # TS name
⋮----
@dataclass
class CoverageEffect
⋮----
"""Coverage effect specification."""
affected_species: str                      # Species affected by coverage
coverage_species: str                      # Species causing coverage
intervals: List[float]                     # Coverage intervals [0, 0.5, 1.0]
slopes: List[float]                        # Slopes for each interval (kcal/mol)
⋮----
# Main Predictor Class
⋮----
class pMuTTPredictor
⋮----
"""
    Wrapper class for pMuTT thermochemistry calculations.

    Simplifies the process of:
    1. Converting DFT data to NASA polynomials
    2. Applying coverage effects
    3. Calculating reaction thermodynamics
    4. Exporting to Chemkin/Cantera formats

    Example usage:
    -------------
    ```python
    from ase.build import molecule

    predictor = pMuTTPredictor()

    # Add gas species from DFT data
    predictor.add_species(
        name='H2O',
        atoms=molecule('H2O'),
        potentialenergy=-14.2209,  # eV
        vib_wavenumbers=[3825.434, 3710.264, 1582.432],
        phase='gas'
    )

    # Generate NASA polynomials
    nasa_species = predictor.generate_nasa_polynomials()

    # Export to Chemkin
    predictor.export_to_chemkin('thermo.dat')

    # Export to Cantera
    predictor.export_to_cantera('mechanism.yaml')
    ```
    """
⋮----
def __init__(self, work_dir: Optional[str] = None, verbose: bool = True)
⋮----
"""
        Initialize pMuTT predictor.

        Args:
            work_dir: Working directory for output files
            verbose: Print detailed information
        """
⋮----
# Storage for species and reactions
⋮----
# Reference species for DFT corrections
⋮----
# Default settings for NASA polynomial generation
self.T_low: float = 200.0    # K
self.T_high: float = 3500.0  # K
self.T_mid: float = 1000.0   # K
⋮----
def _log(self, message: str)
⋮----
"""Print log message if verbose."""
⋮----
# =========================================================================
# Species Management
⋮----
"""
        Add a species to the database.

        Args:
            name: Species name (e.g., 'H2O', 'CO', 'CO*')
            atoms: ASE Atoms object with structure
            potentialenergy: DFT total energy in eV
            vib_wavenumbers: Vibrational frequencies in cm^-1
            symmetrynumber: Symmetry number (1 for no symmetry)
            spin: Electronic spin
            phase: 'gas' or 'surface'
            elements: Element composition {'H': 2, 'O': 1}
            geometry: 'nonlinear', 'linear', or 'atom'
            **kwargs: Additional parameters
        """
# Auto-detect elements from atoms if not provided
⋮----
symbols = atoms.get_chemical_symbols()
elements = dict(Counter(symbols))
⋮----
# Auto-detect geometry
⋮----
geometry = 'atom'
⋮----
geometry = 'linear'
⋮----
species_data = SpeciesData(
⋮----
"""
        Add species from VASP OUTCAR file.

        Args:
            name: Species name
            outcar_path: Path to VASP OUTCAR file
            phase: 'gas' or 'surface'
            **kwargs: Additional parameters
        """
⋮----
atoms = get_atoms(outcar_path)
vib_wavenumbers = get_vib_wavenumbers(outcar_path)
potentialenergy = atoms.get_potential_energy()
⋮----
"""
        Add species from Gaussian log file.

        Args:
            name: Species name
            log_path: Path to Gaussian log file
            phase: 'gas' or 'surface'
            **kwargs: Additional parameters
        """
⋮----
atoms = get_atoms(log_path)
vib_wavenumbers = get_vib_wavenumbers(log_path)
⋮----
def set_reference_species(self, references: List[Dict], T_ref: float = 298.15)
⋮----
"""
        Set reference species for DFT energy corrections.

        Args:
            references: List of reference species dicts
                Example: [
                    {'name': 'H2', 'elements': {'H': 2}, 'HoRT_ref': 0.0},
                    {'name': 'H2O', 'elements': {'H': 2, 'O': 1}, 'HoRT_ref': -96.6}
                ]
            T_ref: Reference temperature (K), default 298.15

        Note:
            For simple usage, you may not need reference species.
            NASA polynomials can be generated directly from StatMech models.
        """
ref_objects = []
⋮----
# Add T_ref if not present
⋮----
ref = Reference(**ref_data)
⋮----
# Coverage Effects
⋮----
"""
        Add coverage effect to a species.

        Args:
            affected_species: Name of species affected by coverage
            coverage_species: Name of species causing coverage
            intervals: Coverage intervals [0, 0.3, 0.7, 1.0]
            slopes: Slopes for each interval in kcal/mol [-8.0, -5.0, -2.0]

        Example:
            # O coverage affects CO adsorption energy
            predictor.add_coverage_effect(
                affected_species='CO*',
                coverage_species='O*',
                intervals=[0.0, 0.3, 0.7, 1.0],
                slopes=[-8.0, -5.0, -2.0]
            )
        """
⋮----
cov_effect = PiecewiseCovEffect(
⋮----
# Add to species data
⋮----
# StatMech Model Generation
⋮----
def generate_statmech_models(self)
⋮----
"""
        Generate statistical mechanics models for all species.

        Returns:
            Dict of StatMech models {name: model}
        """
⋮----
# Choose preset based on phase
⋮----
preset = presets['electronic']
⋮----
preset = presets['idealgas']
else:  # surface
preset = presets['harmonic']
⋮----
# Create StatMech model
statmech_kwargs = {
⋮----
# Add atoms if available
⋮----
# Add vibrations
⋮----
# Create model
model = StatMech(**statmech_kwargs)
⋮----
# Apply coverage effects if any
⋮----
misc_models = [eff['model'] for eff in data.coverage_effects]
⋮----
# NASA Polynomial Generation
⋮----
"""
        Generate NASA polynomials from statistical mechanics models.

        Args:
            T_low: Lower temperature limit (K), default 200
            T_high: Upper temperature limit (K), default 3500
            T_mid: Mid temperature (K), default 1000

        Returns:
            Dict of NASA models {name: model}
        """
⋮----
T_low = T_low or self.T_low
T_high = T_high or self.T_high
T_mid = T_mid or self.T_mid
⋮----
data = self.species_data[name]
⋮----
# Generate NASA polynomial using from_model (replaces deprecated from_statmech)
nasa = Nasa.from_model(
⋮----
# Reaction Management
⋮----
"""
        Add a reaction.

        Args:
            reactants: List of reactant species names
            products: List of product species names
            reactants_stoich: Stoichiometric coefficients (default: all 1)
            products_stoich: Stoichiometric coefficients (default: all 1)
            transition_state: Name of transition state species

        Example:
            predictor.add_reaction(
                reactants=['CO*', 'O*'],
                products=['CO2', '*', '*'],
                transition_state='CO-O*'
            )
        """
# Default stoichiometry
⋮----
reactants_stoich = [1.0] * len(reactants)
⋮----
products_stoich = [1.0] * len(products)
⋮----
# Get species objects
reactant_species = []
⋮----
product_species = []
⋮----
# Transition state
ts_species = None
⋮----
ts_species = self.nasa_models[transition_state]
⋮----
ts_species = self.statmech_models[transition_state]
⋮----
# Create reaction
reaction = Reaction(
⋮----
"""
        Calculate thermodynamic properties for all reactions.

        Args:
            temperature: Temperature in K
            units: Energy units ('kJ/mol', 'eV', 'kcal/mol')

        Returns:
            DataFrame with reaction properties
        """
results = []
⋮----
result = {
⋮----
# Activation energies if TS is available
⋮----
# Export Functions
⋮----
"""
        Export NASA polynomials to Chemkin THERMO.DAT format.

        Args:
            filename: Output filename
            species_list: List of species to export (None = all)
        """
⋮----
# Select species
⋮----
nasa_list = list(self.nasa_models.values())
⋮----
nasa_list = [self.nasa_models[name] for name in species_list]
⋮----
# Export
output_path = os.path.join(self.work_dir, filename)
⋮----
"""
        Export to Cantera YAML format.

        Args:
            filename: Output filename
            gas_species: List of gas species names (None = auto-detect)
            surface_species: List of surface species names (None = auto-detect)
        """
⋮----
# Auto-detect phases if not specified
⋮----
gas_species = [name for name, data in self.species_data.items()
⋮----
surface_species = [name for name, data in self.species_data.items()
⋮----
phases = []
⋮----
# Gas phase
⋮----
gas_nasa = [self.nasa_models[name] for name in gas_species]
gas_phase = IdealGas(
⋮----
# Surface phase
⋮----
surf_nasa = [self.nasa_models[name] for name in surface_species]
surf_phase = StoichSolid(
⋮----
density=2.7e-9  # mol/cm^2, typical for metal surfaces
⋮----
# Convert to CTI string
cti_strings = []
⋮----
cti_str = pmutt_cantera.obj_to_cti(phase)
⋮----
# Write to file
⋮----
"""
        Export thermodynamic data to Excel file.

        Args:
            filename: Output filename
            temperature_range: List of temperatures to evaluate (K)
        """
⋮----
temperature_range = [298.15, 500, 1000, 1500, 2000]
⋮----
# Prepare data
data = []
⋮----
row = {
⋮----
df = pd.DataFrame(data)
⋮----
# Write to Excel
⋮----
# Calculation Methods
⋮----
"""
        Calculate thermodynamic properties for a species at given temperature.

        Args:
            species_name: Name of species
            temperature: Temperature in K
            units: Energy units

        Returns:
            Dict with Cp, H, S, G
        """
⋮----
model = self.nasa_models[species_name]
⋮----
model = self.statmech_models[species_name]
⋮----
def get_species_summary(self) -> pd.DataFrame
⋮----
"""Get summary of all species."""
⋮----
def summary(self) -> str
⋮----
"""Generate a summary of the predictor state."""
lines = ["=" * 60]
⋮----
# Count by phase
gas_count = sum(1 for d in self.species_data.values() if d.phase == 'gas')
surf_count = sum(1 for d in self.species_data.values() if d.phase == 'surface')
⋮----
# Convenience Functions
⋮----
"""
    Quickly generate NASA polynomial from DFT data.

    Args:
        name: Species name
        atoms: ASE Atoms object
        potentialenergy: DFT energy (eV)
        vib_wavenumbers: Vibrational frequencies (cm^-1)
        phase: 'gas' or 'surface'
        **kwargs: Additional parameters

    Returns:
        Nasa object
    """
predictor = pMuTTPredictor(verbose=False)
⋮----
nasa_models = predictor.generate_nasa_polynomials()
⋮----
# Main entry point for testing
⋮----
# Create predictor
predictor = pMuTTPredictor(work_dir='output/pmutt_test')
⋮----
# Test 1: Add gas species
⋮----
# H2
⋮----
# O2
⋮----
# H2O
⋮----
# Test 2: Generate NASA polynomials (no references needed for basic usage)
⋮----
# Calculate properties at 500 K
⋮----
props = predictor.calculate_thermo_at_T(name, 500.0, units='kJ/mol')
⋮----
# Test 3: Add reaction and calculate properties
⋮----
rxn_props = predictor.calculate_reaction_properties(temperature=500.0)
⋮----
# Test 4: Export to different formats
⋮----
# Chemkin
chemkin_file = predictor.export_to_chemkin()
⋮----
# Cantera
cantera_file = predictor.export_to_cantera()
⋮----
# Excel
excel_file = predictor.export_to_excel()
⋮----
# Test 5: Summary
````

## File: core/pathway/__init__.py
````python

````

## File: core/pathway/barrier_predictor.py
````python
"""
Barrier Predictor - 反应能垒预测封装类

这个模块提供了一个完整的封装类，用于使用 Fairchem 的 CatTSunami 框架预测反应能垒。

CatTSunami 使用 Nudged Elastic Band (NEB) 方法来寻找反应过渡态和计算活化能。

使用方式:
    from barrier_predictor import BarrierPredictor
    from fairchem_predictor import FairchemPredictor

    # 初始化预测器（共享 Fairchem 模型）
    fairchem = FairchemPredictor(
        fairchem_root="/path/to/fairchem",
        model_name="uma-s-1p1",
    )

    barrier_predictor = BarrierPredictor(
        fairchem_predictor=fairchem,
    )

    # 方式1: 从反应式预测能垒（需要 CatTSunami 数据库）
    result = barrier_predictor.predict_from_reaction(
        reaction_str="*CH -> *C + *H",
        surface="Pt_111.vasp",
        n_frames=10,
    )
    print(f"Activation energy: {result.activation_energy:.3f} eV")

    # 方式2: 从已知的反应物和产物结构预测能垒
    result = barrier_predictor.predict_from_structures(
        reactant="reactant.vasp",
        product="product.vasp",
        n_frames=10,
    )
    print(f"Activation energy: {result.activation_energy:.3f} eV")

    # 方式3: 手动提供初始 NEB 路径
    result = barrier_predictor.run_neb(
        initial_frames=[atoms1, atoms2, ..., atoms_n],
        fmax=0.1,
    )

作者: Claude
"""
⋮----
# Import LLM-controlled NEB structure preparation
⋮----
@dataclass
class NEBResult
⋮----
"""NEB 计算结果"""
reaction_name: str
reactant: Atoms
product: Atoms
neb_frames: List[Atoms]
n_frames: int
energies: List[float]  # eV
activation_energy_forward: float  # eV (reactant -> product)
activation_energy_reverse: float  # eV (product -> reactant)
reaction_energy: float  # eV (E_product - E_reactant)
transition_state_index: int
transition_state: Atoms
converged: bool
fmax_final: float
optimization_steps: int
trajectory_file: Optional[str] = None
⋮----
def __repr__(self)
⋮----
def summary(self) -> str
⋮----
"""生成结果摘要"""
lines = [
⋮----
@dataclass
class ReactionBarriersResult
⋮----
"""多个反应能垒的结果"""
surface_formula: str
reactions: List[str]
barriers: Dict[str, NEBResult]
output_dir: str
⋮----
# 按正向能垒排序
sorted_reactions = sorted(
⋮----
class BarrierPredictor
⋮----
"""
    反应能垒预测器

    这个类封装了 Fairchem CatTSunami 的功能，提供简洁的 API 用于预测反应能垒。

    Parameters
    ----------
    fairchem_predictor : FairchemPredictor
        Fairchem 预测器实例（用于能量计算）
    fairchem_root : str, optional
        Fairchem 项目根目录（如果不提供 fairchem_predictor）
    work_dir : str, optional
        工作目录
    keep_files : bool, default=True
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """
⋮----
use_llm_controller: bool = False,  # New: enable LLM-controlled NEB preparation
llm_model: str = None,  # New: LLM model for structure preparation
⋮----
# 使用共享的 Fairchem 预测器或创建新的
⋮----
# 延迟导入避免循环依赖
⋮----
# 设置工作目录
⋮----
# 添加 Fairchem 路径
fairchem_src = os.path.join(self.fairchem_root, "src")
⋮----
# 初始化日志
⋮----
# Initialize LLM controller for NEB structure preparation (after logger is set up)
⋮----
# LLM配置（用于智能原子调整）
⋮----
def _init_llm(self)
⋮----
"""初始化LLM用于智能原子调整 - 直接使用OpenAI客户端"""
⋮----
api_key = os.getenv("OPENAI_API_KEY")
base_url = os.getenv("OPENAI_BASE_URL")
model = os.getenv("OPENAI_MODEL", "claude-opus-4-5-20251101")
⋮----
"""
        分析反应物和产物之间的原子差异

        Returns
        -------
        Dict with:
            - 'atoms_to_add': List of (element, count) to add to reactant
            - 'atoms_to_remove': List of (element, count) to remove from reactant
            - 'net_change': Dict of element -> count change
        """
⋮----
reactant_symbols = Counter(reactant.get_chemical_symbols())
product_symbols = Counter(product.get_chemical_symbols())
⋮----
net_change = {}
atoms_to_add = []
atoms_to_remove = []
⋮----
all_elements = set(reactant_symbols.keys()) | set(product_symbols.keys())
⋮----
r_count = reactant_symbols.get(elem, 0)
p_count = product_symbols.get(elem, 0)
diff = p_count - r_count
⋮----
"""
        使用LLM决定要添加的原子的位置 - 提供完整结构信息

        Parameters
        ----------
        atoms : Atoms
            当前结构（反应物）
        adsorbate_indices : List[int]
            吸附分子原子的索引
        atoms_to_add : List[Tuple[str, int]]
            要添加的原子列表 [(element, count), ...]
        reaction_name : str
            反应名称
        reference_atoms : Atoms, optional
            参考结构（产物）
        reference_adsorbate_indices : List[int], optional
            参考结构的吸附分子原子索引

        Returns
        -------
        List[Tuple[str, np.ndarray]]
            要添加的原子列表 [(element, position), ...]
        """
⋮----
# 生成POSCAR格式的结构信息
def atoms_to_poscar_string(structure: Atoms, name: str = "Structure") -> str
⋮----
"""将ASE Atoms对象转换为POSCAR格式字符串"""
lines = [name]
⋮----
# 晶格向量
cell = structure.cell
⋮----
# 元素和数量
symbols = structure.get_chemical_symbols()
unique_elements = []
counts = []
⋮----
# 原子坐标
positions = structure.get_positions()
⋮----
marker = " <-- ADSORBATE" if i in adsorbate_indices else ""
⋮----
# 获取结构信息
reactant_poscar = atoms_to_poscar_string(atoms, "REACTANT")
⋮----
product_poscar = ""
product_adsorbate_info = ""
⋮----
# 为产物生成POSCAR时使用产物的吸附原子索引
def atoms_to_poscar_product(structure: Atoms, ads_indices: List[int], name: str) -> str
⋮----
marker = " <-- ADSORBATE" if i in ads_indices else ""
⋮----
product_poscar = atoms_to_poscar_product(reference_atoms, reference_adsorbate_indices, "PRODUCT")
⋮----
# 提取产物中吸附分子的详细信息
ref_positions = reference_atoms.get_positions()
ref_symbols = reference_atoms.get_chemical_symbols()
product_adsorbate_info = "Product adsorbate atoms (target positions):\n"
⋮----
pos = ref_positions[idx]
⋮----
# 计算表面信息
positions = atoms.get_positions()
symbols = atoms.get_chemical_symbols()
⋮----
# 找到表面原子（非吸附分子的最高z坐标）
surface_atom_z = []
⋮----
surface_z_max = max(surface_atom_z) if surface_atom_z else positions[:, 2].max()
⋮----
# 吸附分子信息
adsorbate_info = "Reactant adsorbate atoms:\n"
⋮----
pos = positions[idx]
⋮----
# 构建详细的LLM提示
prompt = f"""You are a computational chemistry expert helping to set up NEB (Nudged Elastic Band) calculations for surface catalysis.
⋮----
atoms_with_positions = []
⋮----
# 调用LLM
response = self._llm_client.chat.completions.create(
response_text = response.choices[0].message.content
⋮----
# 解析JSON
⋮----
# 尝试提取JSON（处理可能的```json代码块）
cleaned_text = response_text.strip()
⋮----
cleaned_text = re.sub(r'^```(?:json)?\s*', '', cleaned_text)
cleaned_text = re.sub(r'\s*```$', '', cleaned_text)
⋮----
json_match = re.search(r'\[.*\]', cleaned_text, re.DOTALL)
⋮----
result = json.loads(json_match.group())
⋮----
elem = item['element']
pos = np.array(item['position'])
⋮----
# 验证位置合理性
# 使用统一且保守的最小距离阈值，避免对特定元素写死规则
min_dist_threshold = 0.9
⋮----
# 尝试修正位置
corrected_pos = self._correct_atom_position(
⋮----
# 回退：使用参考结构中的位置
⋮----
"""
        验证原子位置是否合理

        Returns
        -------
        Tuple[bool, str]
            (是否有效, 原因说明)
        """
# 检查是否在表面上方
⋮----
# 检查与现有原子的距离
⋮----
distances = np.linalg.norm(positions - pos, axis=1)
min_dist = distances.min()
⋮----
# 检查是否在合理范围内（不要太远）
⋮----
"""
        修正不合理的原子位置

        策略：
        1. 如果有参考结构，使用参考位置但偏移2-3 Å
        2. 否则在吸附物质心上方添加
        """
⋮----
# 尝试从参考结构获取目标位置
target_pos = None
⋮----
target_pos = ref_positions[idx]
⋮----
# 从目标位置向上偏移0.5-0.8 Å（代表预反应状态）
# 关键：偏移量要小，0.5-0.8 Å 足够让NEB找到真正的过渡态
# 太大的偏移（如2.5 Å）会导致NEB找不到正确的反应路径
corrected = target_pos.copy()
# H原子用更小的偏移，因为H-C/H-O键长约1.1 Å
⋮----
corrected[2] += 0.6  # 小偏移，接近产物位置
⋮----
# 回退：在吸附物质心上方
⋮----
ads_center = positions[adsorbate_indices].mean(axis=0)
corrected = ads_center.copy()
# 使用较小的偏移
⋮----
# 最后回退：在原位置上方
corrected = pos.copy()
⋮----
"""
        回退方法：从参考结构中获取原子位置

        如果没有参考结构，则在吸附分子质心附近添加
        """
⋮----
# 计算吸附分子质心
⋮----
# 没有吸附分子索引时，使用表面顶部中心
z_max = positions[:, 2].max()
cell = atoms.cell
ads_center = np.array([cell[0, 0] / 2, cell[1, 1] / 2, z_max + 1.5])
⋮----
new_pos = None
⋮----
# 尝试从参考结构获取位置
⋮----
# 找到参考结构中该元素在吸附分子区域的位置
⋮----
candidate_pos = ref_positions[idx]
# 检查是否已经使用
already_used = any(
⋮----
new_pos = candidate_pos
⋮----
# 如果没找到，使用更智能的定位策略
⋮----
# 策略1：如果有参考结构，找到参考结构中该元素的位置并添加小偏移
⋮----
# 遍历所有参考原子，找到该元素类型
⋮----
candidate_pos = ref_pos.copy()
# 添加小的z偏移（0.5-0.8 Å），使其略高于目标位置
# 这样NEB可以找到真正的过渡态
⋮----
# 策略2：如果还是没找到，在吸附分子上方的合理位置添加
⋮----
# 对于H原子，应该靠近C或O原子（典型键长1.0-1.1 Å）
⋮----
# 找到C或O原子
⋮----
target_atom_pos = positions[idx]
# 在该原子上方1.0 Å处放置H
new_pos = target_atom_pos.copy()
⋮----
# 如果仍然没找到，在质心上方
⋮----
new_pos = ads_center.copy()
⋮----
"""
        使用LLM智能调整反应物和产物结构，使原子数匹配以进行NEB计算

        完全通用的方法，适用于任何涉及原子数变化的反应

        Parameters
        ----------
        reactant : Atoms
            反应物结构
        product : Atoms
            产物结构
        reaction_name : str
            反应名称
        reactant_adsorbate_indices : List[int], optional
            反应物中吸附分子的原子索引（从工作流传递）
        product_adsorbate_indices : List[int], optional
            产物中吸附分子的原子索引（从工作流传递）

        Returns
        -------
        Tuple[Atoms, Atoms, str]
            (调整后的反应物, 调整后的产物, 调整说明)
        """
diff = self._analyze_atom_difference(reactant, product)
⋮----
adjusted_reactant = reactant.copy()
adjusted_product = product.copy()
adjustments = []
⋮----
# 如果产物原子更多（如加氢反应），在反应物中添加原子
⋮----
# 使用LLM确定原子位置
atoms_to_place = self._llm_determine_atom_positions(
⋮----
# 添加原子到结构中
⋮----
new_atom = Atom(elem, position=pos)
⋮----
# 为新添加的原子设置tag（吸附物原子应该是tag=2）
# 获取或创建tags数组
n_atoms = len(adjusted_reactant)
⋮----
# 扩展原有tags
old_tags = list(reactant.get_tags())
new_tags = old_tags + [2] * (n_atoms - len(reactant))  # 新原子tag=2
⋮----
# 创建新tags数组
new_tags = [0] * len(reactant) + [2] * (n_atoms - len(reactant))
⋮----
# 如果反应物原子更多（如脱附反应），在产物中添加原子
⋮----
n_atoms = len(adjusted_product)
⋮----
old_tags = list(product.get_tags())
new_tags = old_tags + [2] * (n_atoms - len(product))  # 新原子tag=2
⋮----
new_tags = [0] * len(product) + [2] * (n_atoms - len(product))
⋮----
adjustment_desc = ", ".join(adjustments) if adjustments else "No adjustment"
⋮----
# 验证调整后原子数匹配
⋮----
# 对齐原子顺序，确保NEB可以正确运行
# NEB要求两个结构的原子顺序完全一致
⋮----
"""
        对齐两个结构的原子顺序，使其适合NEB计算

        通过元素类型和位置匹配原子，确保两个结构的原子索引对应
        """
⋮----
n_atoms = len(reactant)
⋮----
r_symbols = reactant.get_chemical_symbols()
p_symbols = product.get_chemical_symbols()
r_pos = reactant.get_positions()
p_pos = product.get_positions()
⋮----
# 构建成本矩阵：只有相同元素的原子才能匹配
cost_matrix = np.full((n_atoms, n_atoms), 1e10)
⋮----
# 成本 = 位置距离
dist = np.linalg.norm(r_pos[i] - p_pos[j])
⋮----
# 使用匈牙利算法找到最优匹配
⋮----
# 重排产物结构以匹配反应物顺序
new_product = product[col_ind]
⋮----
# 重新计算tags（重要：fairchem需要正确的tags）
# 不能简单重排tags，因为匹配是基于位置而非原子类型
# 需要根据原子类型和z坐标重新确定tags
⋮----
# 从reactant获取表面元素信息
reactant_tags = reactant.get_tags()
reactant_symbols = reactant.get_chemical_symbols()
⋮----
# 识别表面元素（tag=0或1的元素）
surface_elements = set()
⋮----
# 为new_product设置tags
new_tags = []
new_symbols = new_product.get_chemical_symbols()
new_positions = new_product.positions
⋮----
# 找到所有表面原子和吸附原子
surface_atom_indices = []
adsorbate_atom_indices = []
⋮----
# 对表面原子：根据z坐标分bulk和surface
⋮----
surface_z = new_positions[surface_atom_indices, 2]
z_threshold = np.percentile(surface_z, 70)  # 顶部30%为surface
⋮----
new_tags.append(2)  # adsorbate
⋮----
new_tags.append(1)  # surface layer
⋮----
new_tags.append(0)  # bulk layer
⋮----
# 全部是吸附原子（不太可能）
new_tags = [2] * len(new_product)
⋮----
# 保留constraints
⋮----
"""Constrained relaxation of adsorbate + staged atoms to find local minima.

        All adsorbate atoms (including staged) are free; only slab atoms are
        fixed. This ensures both the adsorbate core and staged atoms reach a
        true energy minimum, which is critical for meaningful NEB barriers.

        Modifies reactant and product IN PLACE.
        """
⋮----
# Determine free atoms: all adsorbate indices (includes staged)
free_set = set(ads_indices or []) | set(staged or [])
⋮----
# Use gentler relaxation for endpoints with staged atoms to avoid
# the staged atom collapsing into the bonded position (which would
# eliminate the reactant-product difference needed for NEB).
has_staged = bool(staged)
effective_fmax = max(fmax, 0.3) if has_staged else fmax
effective_steps = min(max_steps, 40) if has_staged else max_steps
⋮----
fixed = [i for i in range(len(atoms)) if i not in free_set]
⋮----
e_before = atoms.get_potential_energy()
⋮----
e_before = None
⋮----
opt = FIRE(atoms, logfile=None)
⋮----
n_steps = opt.get_number_of_steps()
⋮----
e_after = atoms.get_potential_energy()
⋮----
e_after = None
⋮----
"""CatTSunami-style interpolation: linear + iterative overlap correction.

        Based on fairchem/applications/cattsunami/core/autoframe.py::interpolate().
        Uses linear interpolation as initial guess, then iteratively pushes
        apart atoms that are closer than their target distances.

        Args:
            initial_frame: reactant Atoms
            final_frame: product Atoms
            num_frames: total number of frames (including endpoints)
            n_iterations: number of overlap correction iterations
            rate: correction rate per iteration

        Returns:
            list of Atoms frames
        """
⋮----
start_pos = torch.from_numpy(initial_frame.get_positions()).double()
end_pos = torch.from_numpy(final_frame.get_positions()).double()
num_atoms = len(start_pos)
⋮----
# Target distance matrices (linearly interpolated)
start_dist = torch.from_numpy(
end_dist = torch.from_numpy(
⋮----
alpha = torch.linspace(0, 1, num_frames, dtype=torch.float64)
⋮----
# Linear interpolation: frame_i = initial * (1-alpha_i) + final * alpha_i
frames = start_pos.unsqueeze(0) * (1.0 - alpha.view(-1, 1, 1)) + \
⋮----
target_dist = start_dist.unsqueeze(0) * (1.0 - alpha.view(-1, 1)) + \
⋮----
# Build initial Atoms list
def frames_to_atoms(frames_tensor)
⋮----
atoms_list = []
⋮----
a = initial_frame.copy()
⋮----
# Iterative overlap correction
⋮----
atoms_list = frames_to_atoms(frames)
⋮----
frame_dist = []
frame_vec = []
⋮----
d = a.get_all_distances(mic=True)
v = a.get_all_distances(mic=True, vector=True)
⋮----
frame_dist = torch.from_numpy(np.array(frame_dist)).double()
frame_vec = torch.from_numpy(np.array(frame_vec)).double()
⋮----
# Normalize direction vectors
norms = torch.norm(frame_vec, dim=3, keepdim=True) + 1e-8
frame_vec = frame_vec / norms
⋮----
# Only correct atoms that are too close (delta < 0)
delta = torch.clamp(frame_dist - target_dist, max=0)
weight = torch.exp(-(target_dist ** 2) / (0.5 * 2.0 * 2.0))
weight = weight + torch.exp(-(frame_dist ** 2) / (0.5 * 2.0 * 2.0))
delta = delta * weight
⋮----
delta_3d = rate * frame_vec * delta.view(num_frames, num_atoms, num_atoms, 1)
delta_3d = torch.sum(delta_3d, dim=1)  # sum over atom pairs
⋮----
# Freeze endpoints
⋮----
frames = frames - delta_3d
⋮----
# Smooth adjacent frames to avoid large jumps
⋮----
mean_pos = (frames[:-2] + frames[2:]) / 2.0
⋮----
result = frames_to_atoms(frames)
⋮----
# Log min distances for each frame
⋮----
dists = a.get_all_distances(mic=True)
⋮----
min_d = dists.min()
⋮----
def _build_neb_image_calculators(self, n_images: int) -> tuple[list[object], bool]
⋮----
"""Create one calculator per NEB image when the backend supports it.

        ASE's NEB implementation assumes each image owns its own calculator
        instance. Sharing a single FAIRChemCalculator across all images means
        mutable calculator state (`self.atoms`, `self.results`) is reused across
        different endpoint/image objects, which is explicitly discouraged by ASE
        and has shown up as nondeterministic native crashes during GPU NEB runs.

        To keep GPU memory usage reasonable, we do not duplicate the underlying
        ML model. Instead we create lightweight calculator wrappers that all
        point at the same loaded predict unit.
        """
⋮----
base_calc = self.fairchem._calculator
⋮----
calc_cls = type(base_calc)
shared_predictor = getattr(base_calc, "predictor", None)
task_name = getattr(base_calc, "task_name", None)
⋮----
calculators: list[object] = []
⋮----
"""
        运行 NEB 计算

        Uses DyNEB + FIRE with two-phase approach:
          Phase 1: Plain NEB (no climbing) to fmax + delta_fmax_climb (0.2 eV/Å)
          Phase 2: CI-NEB (climbing image) to fmax (0.1 eV/Å)

        Parameters
        ----------
        initial_frames : list of Atoms
            初始 NEB 路径（包括反应物和产物）
        fmax : float, default=0.1
            CI-NEB 力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        delta_fmax_climb: float, default=0.1
            Phase 1 的 fmax 增量，Phase 1 目标为 fmax + delta_fmax_climb
        k : float, default=1.0
            弹簧常数 (eV/Å²)
        climb : bool, default=True
            是否使用爬升像方法
        use_two_phase : bool, default=True
            是否使用两阶段优化（Phase 1无climbing + Phase 2有climbing）
            True: 先以 fmax+delta 稳定路径，再启用 climbing 优化到 fmax（推荐）
            False: 直接使用 CI-NEB 优化到 fmax
        trajectory_file : str, optional
            轨迹文件路径
        reaction_name : str, default="reaction"
            反应名称

        Returns
        -------
        NEBResult
            NEB 计算结果
        """
⋮----
# 复制帧并为每个 image 分配计算器
frames = [atoms.copy() for atoms in initial_frames]
⋮----
# 使用 DyNEB（动态 NEB）
neb = DyNEB(
⋮----
climb=climb if not use_two_phase else False,  # 如果不用两阶段，直接启用climbing
parallel=False,           # 共享 GPU predictor 不支持并行（会数据竞争）
dynamic_relaxation=True,  # 启用动态松弛
scale_fmax=0.0,           # 禁用 fmax 缩放
⋮----
# 设置轨迹文件
⋮----
trajectory_file = os.path.join(
⋮----
# 优化器（FIRE 比 BFGS 对 NEB 更稳定，不易振荡）
optimizer = FIRE(neb, trajectory=trajectory_file)
⋮----
def _run_and_check_convergence(target_fmax: float, step_limit: int, phase_name: str) -> bool
⋮----
step_limit = max(0, int(step_limit))
⋮----
monitor_interval = 2
min_steps_for_stagnation = min(500, max(100, step_limit // 3))
stagnant_checks_limit = 50
improvement_tol = max(0.01, 0.15 * float(target_fmax))
severe_multiplier = 10.0
⋮----
monitor_state = {
⋮----
def _stagnation_monitor() -> None
⋮----
forces = np.array(neb.get_forces(), dtype=float)
⋮----
current_fmax = float(np.max(np.linalg.norm(forces, axis=1)))
⋮----
best_fmax = float(monitor_state["best_fmax"])
⋮----
current_steps = int(optimizer.get_number_of_steps())
⋮----
run_result = optimizer.run(fmax=target_fmax, steps=step_limit)
⋮----
msg = str(exc)
⋮----
# ASE has no detach() method — manually remove from observers list
⋮----
# Fallback for optimizer implementations that don't return bool
⋮----
# 两阶段优化：先不用climbing，再用climbing
# Phase 1 gets at most 60% of total budget; Phase 2 always gets at least 40%
# This ensures the climbing image (critical for saddle point) always runs
phase1_budget = int(max_steps * 0.6)
min_phase2_budget = max_steps - phase1_budget  # at least 40%
⋮----
# 第一阶段：不使用爬升像
⋮----
converged_phase1 = _run_and_check_convergence(
⋮----
converged_phase1 = False
⋮----
optimization_steps = optimizer.get_number_of_steps()
⋮----
# 第二阶段：使用爬升像
# Phase 2 ALWAYS runs — climbing image is essential for locating saddle points
converged = converged_phase1
remaining_steps = max(max_steps - optimization_steps, min_phase2_budget)
⋮----
converged_phase2 = _run_and_check_convergence(
converged = converged_phase2
⋮----
converged = False
⋮----
# 单阶段CI-NEB优化（推荐）
⋮----
converged = _run_and_check_convergence(
⋮----
# 提取结果
energies = []
⋮----
e = frame.get_potential_energy()
⋮----
energies = np.array(energies)
⋮----
# 打印NEB路径能量用于诊断
⋮----
# 找到过渡态 — 仅在 interior images 中寻找（排除端点）
# 与 ASE 的 NEBState.imax 一致: 1 + argmax(energies[1:-1])
# 端点是固定的初末态，不应作为过渡态候选
reactant_energy = energies[0]
product_energy = energies[-1]
⋮----
ts_index = 1 + np.argmax(energies[1:-1])
⋮----
ts_index = np.argmax(energies)
ts_energy = energies[ts_index]
⋮----
# 计算能垒和反应能
activation_energy_forward = max(0.0, ts_energy - reactant_energy)
activation_energy_reverse = max(0.0, ts_energy - product_energy)
reaction_energy = product_energy - reactant_energy
⋮----
# 警告：端点能量高于所有 interior images → 端点可能不是局部极小
⋮----
# 计算最终 fmax
forces = frames[ts_index].get_forces()
fmax_final = np.max(np.linalg.norm(forces, axis=1))
⋮----
# Release CUDA memory held by per-frame calculators.
# NEBResult keeps the Atoms objects but we clear their calc references
# (energies/forces are already extracted above).
⋮----
"""
        从反应物和产物结构预测能垒

        Parameters
        ----------
        reactant : str or Atoms
            反应物结构
        product : str or Atoms
            产物结构
        n_frames : int, default=10
            NEB 帧数
        fmax : float, default=0.1
            CI-NEB 力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        relax_endpoints : bool, default=True
            是否先优化反应物和产物
        reaction_name : str, optional
            反应名称
        reactant_adsorbate_indices : List[int], optional
            反应物中吸附分子的原子索引（从工作流传递，用于LLM智能原子调整）
        product_adsorbate_indices : List[int], optional
            产物中吸附分子的原子索引（从工作流传递，用于LLM智能原子调整）
        reactant_staged_indices : List[int], optional
            反应物中暂驻原子的索引（仅这些原子在端点松弛时可移动）
        product_staged_indices : List[int], optional
            产物中暂驻原子的索引（仅这些原子在端点松弛时可移动）
        previous_step_product : Atoms, optional
            前一步反应的产物结构（用于确保路径连续性）
        previous_step_product_adsorbate_indices : List[int], optional
            前一步产物的吸附质原子索引
        step_index : int, optional
            当前步骤在反应路径中的索引
        use_two_phase : bool, default=True
            是否使用两阶段优化（Phase 1 fmax=0.2 → Phase 2 CI-NEB fmax=0.1）
        delta_fmax_climb : float, default=0.1
            Phase 1 的 fmax 增量
        k : float, default=1.0
            弹簧常数 (eV/Å²)

        Returns
        -------
        NEBResult
            NEB 计算结果
        """
# 读取结构
⋮----
reactant_atoms = read(reactant)
⋮----
reaction_name = os.path.basename(reactant).replace(".vasp", "")
⋮----
reactant_atoms = reactant.copy()
⋮----
reaction_name = "reaction"
⋮----
product_atoms = read(product)
⋮----
product_atoms = product.copy()
⋮----
# 保存原始结构用于能量校正
reactant_atoms_original = reactant_atoms.copy()
product_atoms_original = product_atoms.copy()
adjusted_reactant = None
adjusted_product = None
reactant_result = None
product_result = None
⋮----
# 优化端点
⋮----
r_staged = set(reactant_staged_indices or [])
p_staged = set(product_staged_indices or [])
has_staged = bool(r_staged or p_staged)
⋮----
# Symmetric constrained relaxation of BOTH endpoints.
# Fix: adsorbate core + bottom slab layers
# Free: top 2 slab layers + staged atoms
# This ensures both endpoints are relaxed consistently so
# NEB sees comparable energy landscapes at both ends.
⋮----
relax_fmax = 0.1    # same as NEB fmax
relax_steps = 100   # reasonable budget
⋮----
def _get_slab_free_indices(atoms, staged, adsorbate_idx)
⋮----
"""Only free staged atoms; fix everything else (adsorbate + slab).

                    Conservative approach: staged atoms are the only ones that need
                    to find their local energy minimum. Keeping the slab perfectly
                    fixed ensures both endpoints have identical surface geometry.
                    """
⋮----
def _relax_endpoint(atoms, staged, ads_idx, label)
⋮----
free = _get_slab_free_indices(atoms, staged, ads_idx)
fixed = [i for i in range(len(atoms)) if i not in free]
n_staged = len(staged)
n_surface = len(free) - n_staged
⋮----
opt = LBFGS(atoms, logfile=None)
⋮----
# Relax BOTH endpoints symmetrically
⋮----
# Full unconstrained relaxation (original behavior)
⋮----
reactant_result = self.fairchem.predict_energy(
reactant_atoms = reactant_result.final_structure
⋮----
product_result = self.fairchem.predict_energy(
product_atoms = product_result.final_structure
⋮----
# 检查原子数是否匹配
⋮----
# Use LLM controller if enabled, otherwise use legacy method
⋮----
# Use LLM controller to prepare structures
⋮----
# 保存LLM调整后的初末态结构（用于可视化）
⋮----
adjusted_reactant_file = os.path.join(
adjusted_product_file = os.path.join(
⋮----
# ⭐ FIXED: apply_plan already handles atom reordering
# No need to call align_structures_for_neb which uses Hungarian algorithm
# that would undo the element-type based reordering from apply_plan
⋮----
adjustment_desc = f"LLM-controlled: {plan.reaction_type}"
⋮----
# Fallback to legacy method
⋮----
# Legacy method: rule-based atom adjustment
⋮----
# 对调整后的结构进行优化
⋮----
# 调试：打印tags
⋮----
adj_reactant_result = self.fairchem.predict_energy(
reactant_atoms = adj_reactant_result.final_structure
⋮----
adj_product_result = self.fairchem.predict_energy(
product_atoms = adj_product_result.final_structure
⋮----
# 调整失败，回退到仅计算反应能的模式
⋮----
# 计算反应能但不计算能垒
⋮----
reactant_energy = reactant_result.energy
product_energy = product_result.energy
⋮----
reactant_pred = self.fairchem.predict_energy(reactant_atoms, relax=False)
product_pred = self.fairchem.predict_energy(product_atoms, relax=False)
reactant_energy = reactant_pred.energy
product_energy = product_pred.energy
⋮----
# NOTE: Pre-relaxation of endpoints is intentionally DISABLED.
# Testing showed that pre-relaxing adsorbate/staged atoms causes them
# to collapse into bonded positions, eliminating the energy difference
# between endpoints that NEB needs to find transition states.
# Without pre-relax, NEB finds more non-zero barriers (5/8 vs 2/8)
# with more physically reasonable values (0.07-0.72 eV vs 0-2.1 eV).
⋮----
prepared_pair = prepare_endpoint_pair_for_interpolation(
reactant_atoms = prepared_pair.reactant
product_atoms = prepared_pair.product
⋮----
# 生成插值帧（使用 CatTSunami 风格的线性插值 + 迭代重叠校正）
frames = self._interpolate_with_overlap_correction(
⋮----
# 可视化插值路径
⋮----
# 创建可视化管理器（使用work_dir作为输出目录）
viz_manager = CatalysisVisualizationManager(
⋮----
renderer='tachyon',  # 使用Tachyon高质量渲染
⋮----
interpolated_gif = viz_manager.visualize_neb_trajectory(
⋮----
is_optimized=False,  # 插值路径，非优化后的
⋮----
# 运行 NEB
neb_result = self.run_neb(
⋮----
# 能量校正：处理虚拟原子导致的能量偏移
# 当使用LLM控制器添加虚拟原子时，NEB的反应物能量会改变
# 需要应用能量偏移以使E_act与吸附能量系统一致
⋮----
n_virtual_reactant = len(adjusted_reactant) - len(reactant_atoms_original)
n_virtual_product = (
⋮----
# 虚拟原子添加到反应物
# 从NEB能量减去虚拟原子的贡献（大约每个H原子-2到-4eV）
# 更精确的做法：使用原始反应物能量作为参考
⋮----
# 获取原始反应物的DFT能量（从relax_endpoints的结果）
⋮----
original_reactant_energy = reactant_result.energy
⋮----
original_reactant_energy = self.fairchem.predict_energy(reactant_atoms_original, relax=False).energy
⋮----
# NEB Frame 0 的虚拟原子能量贡献
virtual_energy_contribution = neb_result.energies[0] - (original_reactant_energy + (product_result.energy if hasattr(product_result, 'energy') else 0))
⋮----
# 实际的修正：使用原始反应物能量替代Frame 0
# 然后重新计算能垒
⋮----
# 从NEB能量中移除虚拟原子的贡献
# Frame 0应该对应原始反应物的能量
energies_corrected = np.array(neb_result.energies, dtype=float)
⋮----
# 简单的校正：保持所有帧相对关系，但重新基准化Frame 0
energy_shift = original_reactant_energy - energies_corrected[0]
energies_corrected = energies_corrected + energy_shift
⋮----
# 重新计算能垒
ts_index = np.argmax(energies_corrected)
ts_energy = energies_corrected[ts_index]
⋮----
activation_energy_forward = ts_energy - energies_corrected[0]
activation_energy_reverse = ts_energy - energies_corrected[-1]
reaction_energy = energies_corrected[-1] - energies_corrected[0]
⋮----
# 返回校正后的结果
⋮----
"""
        从反应式预测能垒（使用 CatTSunami AutoFrame）

        Parameters
        ----------
        reaction_str : str
            反应式，如 "*CH -> *C + *H"
        surface : str or Atoms
            表面结构
        n_frames : int, default=10
            NEB 帧数
        n_pdt1_sites : int, default=4
            产物1的位点数
        n_pdt2_sites : int, default=4
            产物2的位点数
        fmax : float, default=0.1
            力收敛标准
        max_steps : int, default=300
            最大优化步数

        Returns
        -------
        NEBResult
            NEB 计算结果
        """
⋮----
# 创建反应对象
reaction = Reaction(
⋮----
# 读取表面
⋮----
surface_atoms = read(surface)
⋮----
surface_atoms = surface.copy()
⋮----
# 创建 Slab 对象
slab = Slab(atoms=surface_atoms)
⋮----
# 1. 生成并优化反应物配置
⋮----
reactant = Adsorbate(
⋮----
reactant_configs = AdsorbateSlabConfig(
⋮----
# 优化反应物
reactant_energies = []
⋮----
result = self.fairchem.predict_energy(
⋮----
best_reactant_idx = np.argmin(reactant_energies)
reactant_system = reactant_configs[best_reactant_idx]
⋮----
# 2. 生成并优化产物配置
⋮----
product1 = Adsorbate(
product2 = Adsorbate(
⋮----
product1_configs = AdsorbateSlabConfig(
⋮----
product2_configs = AdsorbateSlabConfig(
⋮----
# 优化产物
product1_energies = []
⋮----
product2_energies = []
⋮----
# 3. 使用 AutoFrame 生成 NEB 帧
⋮----
af = AutoFrameDissociation(
⋮----
# 4. 运行 NEB（使用第一组帧）
⋮----
"""
        预测多个反应的能垒

        Parameters
        ----------
        surface : str or Atoms
            表面结构
        reactions : list of tuple
            反应列表，每个元组为 (reactant, product, name)
        n_frames : int, default=10
            NEB 帧数
        fmax : float, default=0.1
            力收敛标准
        max_steps : int, default=300
            最大优化步数
        output_dir : str, optional
            输出目录

        Returns
        -------
        ReactionBarriersResult
            多个反应能垒的结果
        """
⋮----
output_dir = os.path.join(self.work_dir, "barriers")
⋮----
surface_formula = surface_atoms.get_chemical_formula()
⋮----
# 预测每个反应
barriers = {}
reaction_names = []
⋮----
# 计算吸附分子的原子索引
surface_atom_count = len(surface_atoms)
reactant_atoms = read(reactant) if isinstance(reactant, str) else reactant
product_atoms = read(product) if isinstance(product, str) else product
reactant_adsorbate_indices = list(range(surface_atom_count, len(reactant_atoms)))
product_adsorbate_indices = list(range(surface_atom_count, len(product_atoms)))
⋮----
result = self.predict_from_structures(
⋮----
# 保存轨迹
traj_file = os.path.join(output_dir, f"{name}_neb.traj")
⋮----
def __del__(self)
⋮----
"""清理临时文件"""
````

## File: core/pathway/candidate_generators.py
````python
"""Candidate generation: two independent exploration modules.

Module 1 — **AgentGuidedPathwayGenerator**
    LLM (or user) supplies a recommended pathway as an ordered list of
    intermediate labels.  This module validates the sequence, fills in
    missing steps, and produces candidate steps to evaluate.
    Fast: only evaluates species on the suggested path.
    Works for ANY reaction on ANY surface — the LLM proposes, UMA evaluates.

Module 2 — **SystematicCRNExplorer**
    Builds a reaction network from scratch using RDKit recursive bond
    operations (same principle as CARE, but element-agnostic).
    Thorough: enumerates all reachable states up to search budget.
    Discovers pathways the LLM might miss.

Both modules produce the same output format: lists of candidate dicts
with ``product_label``, ``product_elements``, ``operation_type``, etc.
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
# ---------------------------------------------------------------------------
# RDKit availability
⋮----
_RDKIT_AVAILABLE: Optional[bool] = None
⋮----
def _check_rdkit() -> bool
⋮----
from rdkit import Chem  # noqa: F401
_RDKIT_AVAILABLE = True
⋮----
_RDKIT_AVAILABLE = False
⋮----
# ═══════════════════════════════════════════════════════════════════════════
# Module 1: Agent-Guided Pathway Generator
⋮----
class AgentGuidedPathwayGenerator
⋮----
"""Generate candidate pathways from LLM / user recommended sequences.

    Input: an ordered list of intermediate labels, e.g.
        ["*CO2", "*COOH", "*CO", "*CHO", "*CH2O", "*CH3O", "*CH3OH", "*CH3", "*CH4"]

    Output: a list of (reactant_label, product_label, operation_type) steps
    that the search engine can directly evaluate with UMA.

    This works for ANY catalytic reaction — the chemical knowledge comes
    from the LLM or user, not from hard-coded rules.
    """
⋮----
"""Convert an ordered intermediate sequence into candidate steps."""
⋮----
steps: List[Dict[str, Any]] = []
⋮----
reactant = intermediates[i]
product = intermediates[i + 1]
r_elements = parse_species_elements(reactant)
p_elements = parse_species_elements(product)
op_type = _infer_operation_type(r_elements, p_elements, reactant, product)
⋮----
"""Generate steps for multiple recommended pathways (competing routes)."""
⋮----
"""Infer elementary step type from element changes."""
r_total = sum(r_elem.values())
p_total = sum(p_elem.values())
h_diff = p_elem.get("H", 0) - r_elem.get("H", 0)
o_diff = p_elem.get("O", 0) - r_elem.get("O", 0)
c_diff = p_elem.get("C", 0) - r_elem.get("C", 0)
⋮----
# Module 2: Systematic CRN Explorer
⋮----
class SystematicCRNExplorer
⋮----
"""Build reaction network from scratch via recursive bond operations.

    Same principle as CARE's dissociate() pipeline, but:
    - Element-agnostic (works for any elements RDKit can represent)
    - Generates bond-breaking, hydrogenation, AND C-C coupling candidates
    - Supports C2–C6 products via coupling of surface intermediates
    - Integrates with CARE bridge when in-domain
    - Produces candidates one layer at a time for the search engine
    """
⋮----
def __init__(self, care_bridge: Any = None, max_carbon: int = 6, max_heavy_atoms: int = 12)
⋮----
self.max_heavy_atoms = max_heavy_atoms  # max non-H atoms in association products
⋮----
# Registry of ALL discovered surface species for association reactions
self._species_registry: Dict[str, Dict[str, Any]] = {}  # label → {elements, smiles}
# Registry of known formulas for isomerization detection
self._formula_registry: Dict[str, List[str]] = {}  # formula → [label1, label2, ...]
⋮----
"""Generate all candidate next-states from a given adsorbate.

        Operations (all element-agnostic):
        1. Bond dissociation (*AB → *A + *B)
        2. Hydrogenation (*A + *H → *AH)
        3. Surface association (*A + *B → *AB, covers: add-O, add-OH, C-C coupling,
           O-O coupling, N-H coupling, water formation, CO₂ formation, etc.)
        4. Eley-Rideal (*A + B(g) → *AB)
        5. Desorption (*A → A(g))
        6. Dissociative adsorption (A₂(g) → 2*A)
        7. Associative desorption (2*A → A₂(g))
        8. PCET (electrochemical: *A + H⁺+e⁻ → *AH, *A + H₂O+e⁻ → *AOH+H⁺)
        9. Isomerization (*AXH ↔ *AHX, 1,2-H migration)
        """
⋮----
elements = parse_species_elements(species_label)
candidates: List[Dict[str, Any]] = []
is_adsorbed = species_label.startswith("*") and "(g)" not in species_label
⋮----
# Gas-phase: molecular adsorption + dissociative adsorption
⋮----
smiles = self._label_to_smiles(species_label)
⋮----
# 1–2. Bond dissociation + Hydrogenation
⋮----
# 3. Surface association (universal: covers add-O, add-OH, C-C, O-O, N-N, ...)
⋮----
# 4. Eley-Rideal
⋮----
# 7. Associative desorption (2*A → A₂(g))
⋮----
# 8. PCET (add OH, add O for electrochem)
⋮----
# 9. Isomerization
⋮----
# Register this species for future association reactions
⋮----
formula_key = _build_label_from_elements(elements, adsorbed=False)
⋮----
"""Generate candidates using RDKit bond operations."""
⋮----
# --- Bond dissociation (surface radical fragments) ---
⋮----
other = frag_labels[1 - i] if len(frag_labels) == 2 else ""
⋮----
# --- Hydrogenation (add H to each unique site) ---
h_products = _enumerate_h_additions(smiles)
⋮----
product_label = _smiles_to_ads_label(product_smi)
⋮----
# Fallback for molecules like CO where RDKit can't directly add H
h_elements = dict(elements)
⋮----
fallback = _lookup_compositions_for_elements(h_elements)
⋮----
fallback = [_build_label_from_elements(h_elements)]
⋮----
# --- PCET ---
hydro_labels = {c["product_label"] for c in candidates if c["operation_type"] == "hydrogenation"}
⋮----
# --- Desorption ---
⋮----
gas_label = label.replace("*", "") + "(g)"
⋮----
def _formula_expand(self, label: str, elements: Dict[str, int]) -> List[Dict[str, Any]]
⋮----
"""Fallback when RDKit is unavailable."""
⋮----
# Hydrogenation
⋮----
# Dehydrogenation
⋮----
dh = dict(elements)
⋮----
# Generic dissociation for each pair of elements
elem_list = [e for e, c in elements.items() if c > 0]
⋮----
frag = dict(elements)
⋮----
# Desorption
⋮----
# --- 3. Universal surface association (*A + *B → *AB) ---
⋮----
"""Generate association candidates by combining this species with ALL
        previously seen surface intermediates.

        This single operation covers:
        - C-C coupling:   *CO + *CO → *C₂O₂
        - Add O:          *CO + *O → *CO₂
        - Add OH:         *CH₃ + *OH → *CH₃OH
        - Water formation: *OH + *H → *H₂O
        - O-O coupling:   *O + *OH → *OOH
        - N-H coupling:   *N + *H → *NH
        - Any other *A + *B combination
        """
⋮----
seen_products: Set[str] = set()
my_heavy = sum(v for k, v in elements.items() if k != "H")
⋮----
other_elem = other_info.get("elements", {})
other_heavy = sum(v for k, v in other_elem.items() if k != "H")
⋮----
# Size guard
total_heavy = my_heavy + other_heavy
⋮----
total_c = elements.get("C", 0) + other_elem.get("C", 0)
⋮----
# Merge elements
product_elem: Dict[str, int] = {}
⋮----
product_label = _build_label_from_elements(product_elem)
⋮----
# Determine operation type from the "other" species
op_type = "association"
⋮----
op_type = "oxidation"
⋮----
op_type = "hydroxylation"
⋮----
op_type = "hydrogenation"
⋮----
op_type = "coupling"
⋮----
# Self-dimerization: *X + *X → *X₂
dimer_heavy = my_heavy * 2
dimer_c = elements.get("C", 0) * 2
⋮----
dimer_elem = {e: c * 2 for e, c in elements.items()}
dimer_label = _build_label_from_elements(dimer_elem)
⋮----
# --- 4. Eley-Rideal (*A + B(g) → *AB) ---
⋮----
"""Eley-Rideal: gas-phase molecule reacts directly with adsorbate."""
⋮----
# Common gas-phase co-reactants
gas_species = {
⋮----
total_heavy = sum(v for k, v in product_elem.items() if k != "H")
⋮----
# --- 6. Dissociative adsorption + 7. Associative desorption ---
⋮----
"""For gas-phase input: molecular adsorption + dissociative adsorption."""
⋮----
ads_label = "*" + species_label.replace("(g)", "").strip().lstrip("*")
⋮----
# Molecular adsorption: A(g) → *A
⋮----
# Dissociative adsorption: A₂(g) → 2*A (for H2, O2, N2)
_DISSOCIATIVE = {
cleaned = species_label.replace("(g)", "").replace("*", "").strip()
⋮----
"""Associative desorption: 2*A → A₂(g) (for *H, *O, *N)."""
⋮----
_ASSOC = {
⋮----
# --- 8. Extended PCET (electrochem: add OH, add O) ---
⋮----
"""Extended PCET: *A + H₂O → *AOH + H⁺+e⁻, *AOH → *AO + H⁺+e⁻."""
⋮----
# PCET add OH: *A + H₂O + e⁻ → *AOH + H⁺
oh_elem = dict(elements)
⋮----
total_heavy = sum(v for k, v in oh_elem.items() if k != "H")
⋮----
# PCET add O: *A + H₂O → *AO + 2H⁺ + 2e⁻
o_elem = dict(elements)
⋮----
total_heavy2 = sum(v for k, v in o_elem.items() if k != "H")
⋮----
# --- 9. Isomerization (1,2-H migration) ---
⋮----
"""Generate isomerization candidates: same formula, different structure.

        Uses RDKit to enumerate H-migration products (move H from one heavy
        atom to an adjacent one).
        """
⋮----
mol = Chem.MolFromSmiles(smiles)
⋮----
original_formula = rdMolDescriptors.CalcMolFormula(mol)
original_canonical = Chem.MolToSmiles(mol)
seen: Set[str] = {original_canonical}
⋮----
mol_h = Chem.AddHs(mol)
# For each H, try moving it to each adjacent heavy atom's neighbor
⋮----
h_idx = atom.GetIdx()
# Find what the H is bonded to
neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
⋮----
source_idx = neighbors[0]
source_atom = mol_h.GetAtomWithIdx(source_idx)
# Find heavy-atom neighbors of source (potential H destinations)
⋮----
dest_idx = dest.GetIdx()
⋮----
# Try moving H: remove H-source bond, add H-dest bond
emol = Chem.RWMol(Chem.RWMol(mol_h))
⋮----
product = Chem.RemoveHs(emol)
psmi = Chem.MolToSmiles(product)
⋮----
# Verify same formula
pm = Chem.MolFromSmiles(psmi)
⋮----
product_label = _smiles_to_ads_label(psmi)
⋮----
# Also check formula registry for known isomers
⋮----
@staticmethod
    def _try_couple_smiles(smi_a: Optional[str], smi_b: Optional[str]) -> Optional[str]
⋮----
"""Try to create a coupled product SMILES by joining two fragments with a C-C bond."""
⋮----
mol_a = Chem.MolFromSmiles(smi_a)
mol_b = Chem.MolFromSmiles(smi_b)
⋮----
# Combine and add C-C bond between first C atoms
combo = Chem.CombineMols(mol_a, mol_b)
emol = Chem.RWMol(combo)
# Find first C in each fragment
c_a = next((a.GetIdx() for a in mol_a.GetAtoms() if a.GetSymbol() == "C"), None)
c_b_offset = mol_a.GetNumAtoms()
c_b = next((a.GetIdx() + c_b_offset for a in mol_b.GetAtoms() if a.GetSymbol() == "C"), None)
⋮----
def _label_to_smiles(self, label: str) -> Optional[str]
⋮----
"""Convert adsorbate label to SMILES."""
⋮----
cleaned = label.replace("*", "").strip()
result = _try_parse_smiles(cleaned)
⋮----
# Router: dispatches to Module 1 and/or Module 2
⋮----
class CandidateGeneratorRouter
⋮----
"""Routes to agent-guided and/or systematic exploration, merges results."""
⋮----
"""Generate candidate next-steps (used by systematic search engine)."""
all_candidates: List[Dict[str, Any]] = []
⋮----
@staticmethod
    def merge_and_deduplicate(candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]
⋮----
seen: Set[str] = set()
deduped: List[Dict[str, Any]] = []
⋮----
key = c.get("product_label", "")
⋮----
# RDKit helpers (shared by SystematicCRNExplorer)
⋮----
def _try_parse_smiles(cleaned: str) -> Optional[str]
⋮----
"""Try to convert a species name to SMILES."""
⋮----
mol = Chem.MolFromSmiles(cleaned)
⋮----
_KNOWN: Dict[str, str] = {
⋮----
# SMILES → readable adsorbate name mapping
# Maps canonical SMILES to commonly used catalysis names
_SMILES_TO_NAME: Dict[str, str] = {
⋮----
"O=CO": "*HCOOH",     # formic acid (CH2O2)
"[O]C=O": "*COOH",    # carboxyl radical (CHO2)
⋮----
"C=O": "*CH2O",       # formaldehyde
"[CH]=O": "*CHO",     # formyl
"CO": "*CH3OH",       # methanol
⋮----
"C": "*CH4",          # methane
⋮----
def _smiles_to_ads_label(smiles: str) -> str
⋮----
"""Convert SMILES to readable adsorbate label like ``*HCOOH``, ``*CHO``."""
⋮----
# Check known name first
⋮----
mol = Chem.MolFromSmiles(smiles, sanitize=True)
⋮----
canonical = Chem.MolToSmiles(mol)
⋮----
# Unsanitized fallback (radical fragments)
mol = Chem.MolFromSmiles(smiles, sanitize=False)
⋮----
ec: Dict[str, int] = {}
⋮----
s = atom.GetSymbol()
⋮----
def _enumerate_bond_breaks(smiles: str) -> List[Tuple[str, List[str], List[Dict[str, int]], str]]
⋮----
"""Break every non-equivalent bond; return fragment labels + element dicts."""
⋮----
mol = Chem.AddHs(mol)
⋮----
ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
⋮----
ranks = list(range(mol.GetNumAtoms()))
⋮----
seen: Set[Tuple[int, int]] = set()
results = []
⋮----
rp = (min(ranks[a1], ranks[a2]), max(ranks[a1], ranks[a2]))
⋮----
sym1 = mol.GetAtomWithIdx(a1).GetSymbol()
sym2 = mol.GetAtomWithIdx(a2).GetSymbol()
desc = f"break {sym1}-{sym2}"
⋮----
emol = Chem.RWMol(mol)
⋮----
frags = Chem.GetMolFrags(emol, asMols=True, sanitizeFrags=False)
⋮----
frag_labels = []
frag_elem_list = []
⋮----
def _enumerate_h_additions(smiles: str) -> List[Tuple[str, str]]
⋮----
"""Add one H to each unique site; return (product_smiles, description).

    Two strategies:
    1. Direct H addition (works for unsaturated atoms like CH2, NH2).
    2. Bond-order reduction + H addition: for fully saturated molecules like
       O=C=O, reduce one double/triple bond to single then add H on both ends.
       This produces e.g. O=C=O + H → OC(=O)[H] (COOH) and O=C(-H)=O (CHO2).
    """
⋮----
results: List[Tuple[str, str]] = []
⋮----
def _try_add(product_mol, desc_prefix: str)
⋮----
product = Chem.RemoveHs(product_mol)
⋮----
pf = rdMolDescriptors.CalcMolFormula(pm)
⋮----
# Strategy 1: direct H addition
⋮----
idx = atom.GetIdx()
sym = atom.GetSymbol()
emol = Chem.RWMol(Chem.AddHs(mol))
h_idx = emol.AddAtom(Chem.Atom(1))
⋮----
# Strategy 2: reduce bond order + add H (for saturated molecules like CO2)
# For each double/triple bond, reduce it by 1 order and add H to one of the atoms
⋮----
bt = bond.GetBondType()
a1_idx = bond.GetBeginAtomIdx()
a2_idx = bond.GetEndAtomIdx()
a1_sym = mol.GetAtomWithIdx(a1_idx).GetSymbol()
a2_sym = mol.GetAtomWithIdx(a2_idx).GetSymbol()
⋮----
new_bt = Chem.BondType.SINGLE
⋮----
new_bt = Chem.BondType.DOUBLE
⋮----
# Add H to atom 1 (reduce bond, H goes to a1)
⋮----
# Add H to atom 2 (reduce bond, H goes to a2)
emol2 = Chem.RWMol(Chem.AddHs(mol))
⋮----
h_idx2 = emol2.AddAtom(Chem.Atom(1))
⋮----
def _lookup_compositions_for_elements(elements: Dict[str, int]) -> List[str]
⋮----
"""Fallback for RDKit H-addition failures (e.g. CO triple bond)."""
⋮----
key = frozenset((k, v) for k, v in elements.items() if v > 0)
_DB: Dict[FrozenSet[Tuple[str, int]], List[str]] = {
smiles_list = _DB.get(key, [])
⋮----
m = Chem.MolFromSmiles(smi)
⋮----
f = rdMolDescriptors.CalcMolFormula(m)
label = f"*{f}"
⋮----
# Common helpers
⋮----
def _build_label_from_elements(elements: Dict[str, int], adsorbed: bool = True) -> str
⋮----
"""Build ``*CHO`` style label from element counts (Hill order)."""
parts = []
⋮----
count = elements.get(elem, 0)
⋮----
count = elements[elem]
⋮----
formula = "".join(parts)
⋮----
def parse_species_elements(label: str) -> Dict[str, int]
⋮----
"""Parse element counts from ``*CHO``, ``CH3OH(g)``, etc."""
cleaned = label.replace("*", "").replace("(g)", "").replace("(s)", "").replace("(l)", "")
⋮----
cleaned = cleaned.split("+")[0].strip()
elements: Dict[str, int] = {}
⋮----
elem = match.group(1)
count = int(match.group(2) or 1)
⋮----
__all__ = [
````

## File: core/pathway/enhanced_neb_validator.py
````python
"""
增强型 NEB 结构验证与修复系统 (Enhanced NEB Validator)

整合并改进现有的验证功能，提供：
1. 基于 atoms.info 标签的通用验证（无需手动输入）
2. 智能结构修复（自动调整不合理结构）
3. 反应类型感知的验证规则
4. 详细的验证报告和修复建议
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
class ValidationStatus(Enum)
⋮----
"""验证状态"""
PASS = "pass"
WARNING = "warning"
FAIL = "fail"
AUTO_FIXED = "auto_fixed"
⋮----
class ReactionType(Enum)
⋮----
"""反应类型"""
HYDROGENATION = "hydrogenation"  # 加氢
DEHYDROGENATION = "dehydrogenation"  # 脱氢
DISSOCIATION = "dissociation"  # 解离
COUPLING = "coupling"  # 偶联
ISOMERIZATION = "isomerization"  # 异构化
UNKNOWN = "unknown"
⋮----
@dataclass
class GeometryIssue
⋮----
"""几何问题描述"""
issue_type: str  # "too_far", "too_close", "unreasonable_bond", etc.
atom_indices: Tuple[int, ...]  # 相关原子索引
description: str
severity: str  # "critical", "warning", "info"
suggestion: str
auto_fixable: bool = False
⋮----
@dataclass
class ValidationReport
⋮----
"""验证报告"""
status: ValidationStatus
reaction_type: ReactionType
initial_valid: bool
final_valid: bool
issues: List[GeometryIssue] = field(default_factory=list)
fixes_applied: List[str] = field(default_factory=list)
recommendations: List[str] = field(default_factory=list)
metadata: Dict[str, Any] = field(default_factory=dict)
⋮----
def to_dict(self) -> Dict
⋮----
def save(self, path: str)
⋮----
class EnhancedNEBValidator
⋮----
"""
    增强型 NEB 结构验证器
    
    特点：
    1. 完全基于 atoms.info 中的标签，无需手动输入
    2. 自动推断反应类型
    3. 智能修复常见问题
    4. 详细的错误报告和建议
    """
⋮----
# 标准共价键长（Å）
STANDARD_BOND_LENGTHS = {
⋮----
# 原子范德华半径（用于检测重叠）
VDW_RADIUS = {
⋮----
# 反应类型特征
REACTION_SIGNATURES = {
⋮----
"critical_distance": 2.0,  # H 应在 2Å 内
⋮----
"atom_change": {},  # 原子数不变
⋮----
def __init__(self, auto_fix: bool = True, fix_attempts: int = 3)
⋮----
"""
        初始化验证器
        
        Parameters
        ----------
        auto_fix : bool
            是否自动修复问题
        fix_attempts : int
            自动修复的最大尝试次数
        """
⋮----
"""
        验证一对 NEB 初末态结构
        
        Parameters
        ----------
        initial : Atoms
            初态结构（必须包含 atoms.info['adsorbate_indices']）
        final : Atoms
            末态结构（必须包含 atoms.info['adsorbate_indices']）
        expected_reaction : str, optional
            期望的反应类型（如 "hydrogenation"）
        
        Returns
        -------
        ValidationReport
            验证报告，包含问题和修复建议
        """
issues = []
fixes = []
⋮----
# 1. 检查必要的元数据
⋮----
# 2. 推断反应类型
reaction_type = self._infer_reaction_type(initial, final, expected_reaction)
⋮----
# 3. 验证初态结构
initial_issues = self._validate_single_structure(initial, "initial")
⋮----
# 4. 验证末态结构
final_issues = self._validate_single_structure(final, "final")
⋮----
# 5. 验证初末态一致性
consistency_issues = self._validate_consistency(initial, final, reaction_type)
⋮----
# 6. 反应类型特定的验证
specific_issues = self._validate_reaction_specific(initial, final, reaction_type)
⋮----
# 7. 自动修复（如果启用）
⋮----
# 重新验证
⋮----
# 8. 确定最终状态
critical_issues = [i for i in issues if i.severity == "critical"]
warnings = [i for i in issues if i.severity == "warning"]
⋮----
status = ValidationStatus.FAIL
⋮----
status = ValidationStatus.WARNING
⋮----
status = ValidationStatus.AUTO_FIXED
⋮----
status = ValidationStatus.PASS
⋮----
# 9. 生成建议
recommendations = self._generate_recommendations(issues, reaction_type)
⋮----
"""从结构变化推断反应类型"""
⋮----
# 计算原子数变化
initial_symbols = initial.get_chemical_symbols()
final_symbols = final.get_chemical_symbols()
⋮----
initial_counts = Counter(initial_symbols)
final_counts = Counter(final_symbols)
⋮----
# 检测 H 原子变化
h_change = final_counts.get('H', 0) - initial_counts.get('H', 0)
⋮----
# 检测是否解离（吸附物原子数增加）
initial_adsorbate = len(initial.info.get('adsorbate_indices', []))
final_adsorbate = len(final.info.get('adsorbate_indices', []))
⋮----
"""验证单个结构"""
⋮----
positions = atoms.get_positions()
symbols = atoms.get_chemical_symbols()
adsorbate_indices = atoms.info.get('adsorbate_indices', [])
⋮----
# 1. 检查原子重叠
⋮----
dist = np.linalg.norm(positions[idx1] - positions[idx2])
⋮----
min_dist = self.VDW_RADIUS.get(elem1, 1.5) + self.VDW_RADIUS.get(elem2, 1.5)
min_dist *= 0.5  # 允许一定程度的重叠
⋮----
# 2. 检查吸附物-表面距离
surface_indices = [i for i in range(len(atoms)) if i not in adsorbate_indices]
⋮----
surface_z = max(positions[i, 2] for i in surface_indices)
⋮----
height = positions[idx, 2] - surface_z
⋮----
# 3. 检查吸附物内部键长
⋮----
elem_pair = tuple(sorted([symbols[idx1], symbols[idx2]]))
⋮----
"""验证初末态一致性"""
⋮----
# 1. 检查表面原子数是否守恒
initial_surface = len([i for i in range(len(initial))
final_surface = len([i for i in range(len(final))
⋮----
# 2. 检查是否有不合理的原子数变化
⋮----
# 排除表面元素（通常是金属）
surface_elem = max(initial_counts, key=initial_counts.get)
⋮----
change = final_counts.get(elem, 0) - initial_counts.get(elem, 0)
⋮----
"""反应类型特定的验证"""
⋮----
signature = self.REACTION_SIGNATURES.get(reaction_type)
⋮----
# 加氢反应特定检查
⋮----
# 脱氢反应特定检查
⋮----
"""检查加氢反应的几何结构"""
⋮----
# 在末态中，H 应该靠近 C/O/N
final_pos = final.get_positions()
⋮----
final_adsorbate = final.info.get('adsorbate_indices', [])
⋮----
# 找到新增的 H
initial_counts = {}
⋮----
final_h_atoms = [i for i in final_adsorbate if final_symbols[i] == 'H']
⋮----
# 找到最近的 C/O/N
min_dist = float('inf')
min_idx = -1
min_elem = ""
⋮----
dist = np.linalg.norm(final_pos[h_idx] - final_pos[idx])
⋮----
min_dist = dist
min_idx = idx
min_elem = final_symbols[idx]
⋮----
"""检查脱氢反应的几何结构"""
⋮----
# 在初态中，应该有一个合理的 C-H/O-H/N-H 键
initial_pos = initial.get_positions()
⋮----
initial_adsorbate = initial.info.get('adsorbate_indices', [])
⋮----
h_atoms = [i for i in initial_adsorbate if initial_symbols[i] == 'H']
⋮----
dist = np.linalg.norm(initial_pos[h_idx] - initial_pos[idx])
min_dist = min(min_dist, dist)
⋮----
if min_dist > 1.5:  # C-H 键通常 ~1.1 Å
⋮----
"""尝试自动修复问题"""
applied_fixes = []
⋮----
# 修复：将 H 原子移动到目标原子附近
⋮----
initial = fixed_initial
⋮----
def _fix_h_position(self, atoms: Atoms, issue: GeometryIssue) -> Tuple[Atoms, str]
⋮----
"""修复 H 原子位置"""
⋮----
atoms_fixed = atoms.copy()
positions = atoms_fixed.get_positions()
⋮----
target_pos = positions[target_idx]
⋮----
# 将 H 放在目标原子附近，稍微偏上方
new_pos = target_pos + np.array([0.5, 0.5, 1.5])  # 1.5 Å above and offset
⋮----
"""生成修复建议"""
recommendations = []
⋮----
critical = [i for i in issues if i.severity == "critical"]
⋮----
for issue in critical[:3]:  # 只显示前3个
⋮----
# 反应类型特定建议
⋮----
# =============================================================================
# 便捷函数
⋮----
"""
    便捷函数：验证 NEB 初末态
    
    Example:
        >>> from ase.io import read
        >>> initial = read("initial.vasp")
        >>> final = read("final.vasp")
        >>> report = validate_neb_endpoints(initial, final, "hydrogenation")
        >>> print(report.status)
        >>> if report.status != ValidationStatus.PASS:
        ...     print(report.recommendations)
    """
validator = EnhancedNEBValidator(auto_fix=auto_fix)
````

## File: core/pathway/fairchem_predictor.py
````python
"""
Fairchem Predictor - 催化剂自由能预测封装类

这个模块提供了一个完整的封装类，用于使用 Fairchem 模型预测催化剂表面的吸附能和自由能。

Fairchem 提供了多个预训练模型（如 UMA、eSCN、GemNet 等）用于催化剂系统的能量预测。

使用方式:
    from fairchem_predictor import FairchemPredictor

    # 初始化预测器
    predictor = FairchemPredictor(
        fairchem_root="/path/to/fairchem",
        model_name="uma-s-1p1",  # 或 "uma-m-1p1"
        use_gpu=True,
    )

    # 方式1: 预测单个结构的能量
    result = predictor.predict_energy(
        structure="surface_with_adsorbate.vasp",
        relax=True,  # 是否进行结构优化
    )
    print(f"Energy: {result.energy:.3f} eV")

    # 方式2: 预测吸附能
    ads_result = predictor.predict_adsorption_energy(
        surface="clean_surface.vasp",
        adsorbate="*CO",  # 或提供分子结构
        num_sites=10,  # 尝试多少个吸附位点
    )
    print(f"Adsorption energy: {ads_result.adsorption_energy:.3f} eV")

    # 方式3: 批量预测多个吸附质的能量
    pathway_result = predictor.predict_pathway_energies(
        surface="surface.vasp",
        adsorbates=["*O", "*OH", "*OOH"],
        num_sites=5,
    )
    for ads, energy in pathway_result.adsorbate_energies.items():
        print(f"{ads}: {energy:.3f} eV")

作者: Claude
"""
⋮----
@dataclass
class EnergyResult
⋮----
"""单个结构的能量预测结果"""
structure_name: str
initial_structure: Atoms
final_structure: Atoms
energy: float  # eV
forces: np.ndarray  # eV/Å
relaxed: bool
converged: bool
optimization_steps: int
output_file: Optional[str] = None
⋮----
def __repr__(self)
⋮----
@dataclass
class AdsorptionResult
⋮----
"""吸附能预测结果"""
adsorbate: str
surface_energy: float  # eV
adsorbate_surface_energy: float  # eV
adsorbate_reference_energy: float  # eV
adsorption_energy: float  # eV
num_sites_tried: int
best_site_index: int
best_configuration: Atoms
all_configurations: List[Atoms] = field(default_factory=list)
all_energies: List[float] = field(default_factory=list)
is_anomalous: bool = False
⋮----
@dataclass
class PathwayResult
⋮----
"""反应路径能量预测结果"""
surface_formula: str
⋮----
adsorbates: List[str]
adsorbate_energies: Dict[str, float]  # adsorbate -> energy (eV)
adsorption_energies: Dict[str, float]  # adsorbate -> E_ads (eV)
best_configurations: Dict[str, Atoms]
relative_energies: Dict[str, float]  # relative to most stable
output_dir: str
⋮----
def summary(self) -> str
⋮----
"""生成结果摘要"""
lines = [
⋮----
# 按吸附能排序
sorted_ads = sorted(self.adsorption_energies.items(),
⋮----
total_e = self.adsorbate_energies[ads]
rel_e = self.relative_energies[ads]
⋮----
class FairchemPredictor
⋮----
"""
    Fairchem 自由能预测器

    这个类封装了 Fairchem 项目的功能，提供简洁的 API 用于预测催化剂表面的
    吸附能和自由能。

    Parameters
    ----------
    fairchem_root : str
        Fairchem 项目的根目录路径
    model_name : str, default="uma-s-1p1"
        使用的模型名称。可选:
        - "uma-s-1p1": UMA small model (快速)
        - "uma-m-1p1": UMA medium model (更准确)
        - 或提供本地模型文件的完整路径
    model_path : str, optional
        本地模型文件路径。如果提供，将使用本地模型而不是从 Hugging Face 下载
    task_name : str, default="oc20"
        任务类型。催化应用使用 "oc20"
    use_gpu : bool, default=True
        是否使用 GPU
    device : str, optional
        指定设备 ("cuda" 或 "cpu")
    inference_mode : str, default="default"
        推理模式: "default" 或 "turbo" (更快)
    work_dir : str, optional
        工作目录
    keep_files : bool, default=True
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """
⋮----
# 原子参考能量 (eV) - 用于计算吸附能
ATOMIC_REFERENCE_ENERGIES = {
⋮----
env_device = str(os.getenv("CATDT_FAIRCHEM_DEVICE", "")).strip().lower()
requested_device = str(device or ("cuda" if use_gpu else "cpu")).strip().lower()
resolved_device = env_device or requested_device or "cuda"
⋮----
resolved_device = requested_device or "cuda"
⋮----
# 设置工作目录
⋮----
# 添加 Fairchem 路径到系统路径
fairchem_src = os.path.join(self.fairchem_root, "src")
⋮----
# 初始化日志
⋮----
# 延迟加载模型（在第一次使用时）
⋮----
def _load_model(self)
⋮----
"""延迟加载模型"""
⋮----
# PyTorch 2.8+ defaults weights_only=True; UMA checkpoints need `slice` whitelisted
⋮----
# 🔧 确保使用 deps/fairchem 的版本，而不是环境中已安装的旧版本
fairchem_src_path = os.path.join(self.fairchem_root, "src")
fairchem_pkg_path = os.path.join(fairchem_src_path, "fairchem")
⋮----
# 避免在长时间多 episode 运行中反复清空 fairchem 相关模块。
# 这种“硬重置”会触发不可预测的导入态损坏（出现随机 sys/sys、tuple callable 等异常）。
# 仅在显式开启时执行。
hard_reset = str(os.getenv("CATDT_FAIRCHEM_HARD_RESET", "0")).strip().lower() in {
def _reset_fairchem_modules(full: bool = False) -> None
⋮----
local_src_abs = os.path.abspath(fairchem_src_path)
⋮----
mod_file = getattr(mod, "__file__", "") or ""
mod_file_abs = os.path.abspath(mod_file) if mod_file else ""
⋮----
def _import_fairchem_components() -> Tuple[Any, Any, Any]
⋮----
# 检查是否使用本地模型
⋮----
# 使用本地模型文件
⋮----
# 优先使用本地模型，避免网络下载阻塞；失败后再尝试 HuggingFace。
local_checkpoints = [
⋮----
f"deps/fairchem_models/{self.model_name}.pt",  # 本地UMA模型
⋮----
existing_checkpoints = []
⋮----
fallback_model = existing_checkpoints[0] if existing_checkpoints else None
local_errors: List[str] = []
⋮----
max_local_attempts = max(
⋮----
loaded = False
last_err = None
⋮----
loaded = True
⋮----
last_err = err
⋮----
local_error = "; ".join(local_errors) if local_errors else None
⋮----
# 创建计算器
⋮----
def _prepare_slab_for_adsorbate_placement(self, atoms: Atoms) -> Atoms
⋮----
"""
        准备表面用于吸附质放置

        确保表面满足 Fairchem Slab 的要求：
        - Tagged（表面原子 tag=1，体相原子 tag=0）
        - Tiled（xy 方向 >= 8 埃）
        - 添加了 FixAtoms constraints

        Parameters
        ----------
        atoms : Atoms
            原始表面结构

        Returns
        -------
        Atoms
            准备好的表面结构
        """
⋮----
prepared_atoms = atoms.copy()
⋮----
# 1. Tag 原子
# 简单策略：基于 z 坐标区分表面和体相
# 最上面 2 层为表面原子（tag=1），其余为体相原子（tag=0）
z_coords = prepared_atoms.positions[:, 2]
z_max = z_coords.max()
⋮----
# 识别层（通过 z 坐标分组）
z_unique = np.unique(np.round(z_coords, decimals=2))
z_sorted = np.sort(z_unique)[::-1]  # 从上到下排序
⋮----
# 标记前 2-3 层为表面原子
n_surface_layers = min(3, len(z_sorted))
surface_z_threshold = z_sorted[n_surface_layers-1] if len(z_sorted) > 0 else z_max - 2.0
⋮----
tags = np.zeros(len(prepared_atoms), dtype=int)
⋮----
tags[i] = 1  # 表面原子
⋮----
tags[i] = 0  # 体相原子
⋮----
# 确保 PBC 正确设置（表面应该在 xy 方向周期）
⋮----
# 2. 检查是否需要 tile
cell = prepared_atoms.cell
min_ab = 8.0
⋮----
# 需要 tile
nx = int(np.ceil(min_ab / cell[0, 0]))
ny = int(np.ceil(min_ab / cell[1, 1]))
⋮----
prepared_atoms = prepared_atoms * (nx, ny, 1)
# 确保 PBC 仍然正确
⋮----
# 3. 添加 constraints（固定体相原子）
indices_to_fix = [i for i, tag in enumerate(prepared_atoms.get_tags()) if tag == 0]
⋮----
constraint = FixAtoms(indices=indices_to_fix)
⋮----
"""
        获取吸附物的结合原子信息（z坐标最低的原子）

        Parameters
        ----------
        config : Atoms
            表面+吸附物的结构
        n_surface : int
            表面原子数量

        Returns
        -------
        Tuple[int, np.ndarray]
            (结合原子在吸附物中的索引, 结合原子的位置)
        """
adsorbate_indices = list(range(n_surface, len(config)))
⋮----
adsorbate_positions = config.positions[adsorbate_indices]
binding_atom_local_idx = np.argmin(adsorbate_positions[:, 2])
binding_atom_global_idx = adsorbate_indices[binding_atom_local_idx]
binding_atom_pos = config.positions[binding_atom_global_idx].copy()
⋮----
"""
        检查结合原子是否从目标位点漂移

        Parameters
        ----------
        actual_pos : np.ndarray
            实际位置
        target_pos : np.ndarray
            目标位置
        xy_threshold : float
            xy平面内的漂移阈值（埃）

        Returns
        -------
        Tuple[bool, float]
            (是否漂移超出阈值, xy平面内的漂移距离)
        """
xy_drift = np.sqrt((actual_pos[0] - target_pos[0])**2 +
is_drifted = xy_drift > xy_threshold
⋮----
"""
        使用xy方向软约束优化结构

        通过在每几步后重置结合原子的xy位置来实现软约束

        Parameters
        ----------
        config : Atoms
            初始结构
        n_surface : int
            表面原子数量
        target_xy : np.ndarray
            目标xy位置 [x, y]
        fmax : float
            力收敛标准
        max_steps : int
            最大优化步数

        Returns
        -------
        Atoms
            优化后的结构
        """
⋮----
atoms = config.copy()
⋮----
# 找到结合原子索引
adsorbate_indices = list(range(n_surface, len(atoms)))
⋮----
adsorbate_positions = atoms.positions[adsorbate_indices]
⋮----
# 保存原始约束
original_constraints = atoms.constraints.copy() if atoms.constraints else []
⋮----
# 添加xy方向约束（只允许z方向移动）
xy_constraint = FixCartesian(binding_atom_global_idx, mask=[True, True, False])
⋮----
# 先将结合原子移动到目标xy位置
current_pos = atoms.positions[binding_atom_global_idx].copy()
⋮----
# 优化
opt = LBFGS(atoms)
⋮----
# 恢复原始约束
⋮----
"""
        手动放置吸附质并预测吸附能（后备方法）

        当 Fairchem 的 AdsorbateSlabConfig 失败时使用此方法
        或当提供了固定吸附位点时直接使用
        """
⋮----
# 如果提供了固定吸附位点，只使用这个位点
⋮----
sites = [fixed_site_position]
⋮----
# 在表面最高的几个原子周围生成位点
z_coords = surface_atoms.positions[:, 2]
⋮----
# 找到表面原子（z > z_max - 1.5 埃）
surface_atom_indices = np.where(z_coords > z_max - 1.5)[0]
⋮----
# 生成吸附位点：在表面原子上方
sites = []
⋮----
site = surface_atoms.positions[idx].copy()
site[2] = z_max + 2.0  # 在表面上方 2 埃
⋮----
# 如果表面原子数量不足，在它们周围生成更多位点
⋮----
# 在现有位点周围添加偏移位点
⋮----
base_site = sites[i % len(sites)].copy()
# 添加随机偏移
⋮----
# 优化每个配置
all_energies = []
all_configs = []
⋮----
# 创建配置
config = surface_atoms.copy()
ads_copy = adsorbate_atoms.copy()
⋮----
# 将吸附质的结合原子（z坐标最低的原子）移动到位点
# 这确保了不同吸附质使用相同的结合位点
ads_positions = ads_copy.get_positions()
binding_atom_idx = np.argmin(ads_positions[:, 2])  # z坐标最低的原子
binding_atom_pos = ads_positions[binding_atom_idx]
⋮----
# 计算需要的平移量，使结合原子位于site位置
translation = site - binding_atom_pos
⋮----
# 合并（确保 PBC 一致）
n_surface = len(config)
⋮----
# 确保整个系统的 PBC 正确设置
⋮----
# 正确设置tags：保留表面原有的tags，只为吸附质设置tag=2
⋮----
# 使用表面原有的tags
surface_tags = list(surface_atoms.get_tags())
⋮----
# 如果没有tags，从constraints推断
⋮----
fixed_indices = []
⋮----
surface_tags = [0 if i in fixed_indices else 1 for i in range(n_surface)]
⋮----
# 没有constraints和tags，假设底部原子是bulk
surface_tags = [1] * n_surface
⋮----
# 吸附质原子tag=2
tags = surface_tags + [2] * len(ads_copy)
⋮----
# 验证tags设置
set_tags = config.get_tags()
⋮----
result = self.predict_energy(
⋮----
# 验证优化后tags
final_tags = result.final_structure.get_tags()
⋮----
final_config = result.final_structure
final_energy = result.energy
⋮----
# 检查结合原子是否从目标位点漂移（仅在使用固定位点时检查）
⋮----
# 重新优化：使用xy约束
constrained_config = self._optimize_with_xy_constraint(
⋮----
config=config,  # 使用初始配置（吸附物在目标位点）
⋮----
target_xy=site[:2],  # 只用xy坐标
⋮----
# 计算约束优化后的能量
⋮----
constrained_energy = constrained_config.get_potential_energy()
⋮----
# 检查约束后的位置
⋮----
# 使用约束优化的结果（保持位点一致性比能量更重要）
final_config = constrained_config
final_energy = constrained_energy
⋮----
# 恢复tags
⋮----
pass  # tags应该已经在约束优化中保留
⋮----
# 找到最低能量配置
best_idx = np.argmin(all_energies)
best_energy = all_energies[best_idx]
best_config = all_configs[best_idx]
⋮----
# 计算吸附能
adsorption_energy = best_energy - surface_energy - adsorbate_ref_energy
⋮----
"""
        预测单个结构的能量

        Parameters
        ----------
        structure : str or Atoms
            输入结构（文件路径或 ASE Atoms 对象）
        relax : bool, default=True
            是否进行结构优化
        fmax : float, default=0.05
            力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        optimizer : str, default="LBFGS"
            优化器类型 ("LBFGS" 或 "BFGS")
        output_file : str, optional
            输出文件路径

        Returns
        -------
        EnergyResult
            能量预测结果
        """
⋮----
# 读取结构
⋮----
atoms = read(structure)
structure_name = os.path.basename(structure)
⋮----
atoms = structure.copy()
structure_name = atoms.get_chemical_formula()
⋮----
# 确保 PBC 设置一致（FAIRChemCalculator 要求）
# PBC 必须是全 True 或全 False
pbc = atoms.pbc
⋮----
# 如果有部分为 True，则全部设为 True
⋮----
initial_structure = atoms.copy()
⋮----
# 保存原始tags（fairchem计算器可能会修改tags）
original_tags = atoms.get_tags().copy() if atoms.has('tags') else None
⋮----
# 结构优化
converged = False
optimization_steps = 0
⋮----
opt = BFGS(atoms)
⋮----
converged = True
optimization_steps = opt.get_number_of_steps()
⋮----
optimization_steps = max_steps
⋮----
# 检查优化后的tags
⋮----
# 恢复原始tags（必须在获取能量前恢复，因为fairchem可能依赖tags）
⋮----
# 获取最终能量和力
⋮----
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
⋮----
# 如果失败，尝试重新设置calculator
⋮----
# 保存结果
⋮----
"""
        预测吸附能

        Parameters
        ----------
        surface : str or Atoms
            清洁表面结构
        adsorbate : str or Atoms
            吸附质（SMILES 字符串或 Atoms 对象）
        num_sites : int, default=10
            尝试的吸附位点数量
        mode : str, default="random_site_heuristic_placement"
            位点生成模式
        fmax : float, default=0.05
            力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        detect_anomalies : bool, default=True
            是否检测异常（解吸、解离等）
        fixed_site_position : np.ndarray, optional
            固定的吸附位点坐标 (x, y, z)，如果提供则只在此位点添加吸附物
        skip_surface_relaxation : bool, default=False
            是否跳过表面优化（当表面已经优化过时使用）

        Returns
        -------
        AdsorptionResult
            吸附能预测结果
        """
⋮----
# 读取表面结构
⋮----
surface_atoms = read(surface)
⋮----
surface_atoms = surface.copy()
⋮----
# 1. 优化清洁表面（如果需要）
⋮----
# 保留原始surface，包括tags
surface_atoms_with_tags = surface_atoms.copy()
⋮----
# 计算表面能量
surface_atoms_copy = surface_atoms.copy()
⋮----
surface_energy = surface_atoms_copy.get_potential_energy()
⋮----
# 创建result，确保final_structure保留tags
surface_result = EnergyResult(
⋮----
final_structure=surface_atoms_with_tags,  # 保留tags的surface
⋮----
tags = surface_atoms_with_tags.get_tags()
⋮----
surface_result = self.predict_energy(
surface_energy = surface_result.energy
⋮----
# 2. 创建吸附质对象
adsorbate_obj = None
⋮----
adsorbate_obj = Adsorbate(
adsorbate_atoms = adsorbate_obj.atoms
⋮----
adsorbate_atoms = adsorbate.copy()
⋮----
# 3. 计算吸附质参考能量
adsorbate_symbols = adsorbate_atoms.get_chemical_symbols()
adsorbate_ref_energy = sum([
⋮----
# 如果提供了固定吸附位点，直接使用手动放置方法
⋮----
num_sites=1,  # 只用一个位点
⋮----
# 为了使用 AdsorbateSlabConfig，我们需要确保表面满足要求：
# - 被 tagged（表面原子 tag=1，体相原子 tag=0）
# - 被 tiled（xy 方向 >= 8 埃）
# - 添加了 constraints
slab_atoms = self._prepare_slab_for_adsorbate_placement(
⋮----
# 创建 Slab 对象
⋮----
slab = Slab(slab_atoms=slab_atoms)
⋮----
# 如果 Slab 创建失败，使用手动放置方法
⋮----
# 创建吸附质对象（如果还没有）
⋮----
ads_for_config = adsorbate_obj
⋮----
# 需要创建 Adsorbate 对象
⋮----
ads_for_config = Adsorbate(
⋮----
adsorbate_binding_indices=[0]  # 假设第一个原子是绑定原子
⋮----
# 使用手动方法
⋮----
configs = AdsorbateSlabConfig(
⋮----
# 使用手动放置方法
⋮----
# 5. 优化所有配置
⋮----
# 6. 找到最低能量配置
⋮----
# 7. 计算吸附能: E_ads = E(ads+surf) - E(surf) - E(ads_ref)
⋮----
# 8. 异常检测
is_anomalous = False
⋮----
detector = DetectTrajAnomaly(
⋮----
is_anomalous = (
⋮----
"""
        预测反应路径上多个吸附质的能量

        Parameters
        ----------
        surface : str or Atoms
            催化剂表面结构（clean表面，不含吸附物）
        adsorbates : list of str or Atoms
            吸附质列表
        num_sites : int, default=5
            每个吸附质尝试的位点数（如果提供fixed_adsorption_site则忽略）
        fmax : float, default=0.05
            力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        output_dir : str, optional
            输出目录
        fixed_adsorption_site : np.ndarray, optional
            固定的吸附位点坐标 (x, y, z)，如果提供则所有吸附物都使用此位点
        skip_surface_relaxation : bool, default=False
            是否跳过该函数内的 clean-surface relax（已优化表面可启用以提速/稳健）

        Returns
        -------
        PathwayResult
            反应路径能量预测结果
        """
⋮----
output_dir = os.path.join(self.work_dir, "pathway_results")
⋮----
# 读取表面
⋮----
surface_formula = surface_atoms.get_chemical_formula()
⋮----
# 1. 处理清洁表面
⋮----
surface_relaxed = surface_atoms.copy()
surface_calc = surface_atoms.copy()
⋮----
surface_energy = float(surface_calc.get_potential_energy())
⋮----
surface_relaxed = surface_result.final_structure
⋮----
# 保存清洁表面
surface_file = os.path.join(output_dir, "surface_relaxed.vasp")
⋮----
# 2. 预测每个吸附质
adsorbate_energies = {}
adsorption_energies = {}
best_configurations = {}
failed_adsorbates: List[Tuple[str, str]] = []
⋮----
# 跟踪实际的参考吸附位点（从第一个非"*"吸附物的结果中获取）
# 这确保了所有后续吸附物使用相同的位点，即使优化过程中发生了微小漂移
actual_reference_site = None
n_surface_atoms = len(surface_relaxed)
⋮----
# 特殊处理：如果吸附物是 "*"
ads_str = str(ads).strip()
⋮----
ads_name = "*"
⋮----
# 判断是吸附的起始态还是脱附的产物态
⋮----
# 第一个中间体是"*"：吸附反应的起始态（clean surface）
⋮----
ads_file = os.path.join(output_dir, f"{ads_name}_relaxed.vasp")
⋮----
# 中间或末尾的"*"：脱附反应的产物态
# 需要将前一个中间体的分子放在表面上方3-4埃（代表气相/弱吸附态）
⋮----
# 获取前一个中间体
prev_ads_name = str(adsorbates[i-1]).strip()
⋮----
prev_config = best_configurations[prev_ads_name]
⋮----
# 识别前一个配置中的吸附分子原子
n_surface = len(surface_relaxed)
adsorbate_indices = list(range(n_surface, len(prev_config)))
⋮----
# 创建脱附态：分子在表面上方3.5埃
desorbed_config = surface_relaxed.copy()
⋮----
# 提取吸附分子（需要先复制再修改位置）
ads_atoms_list = []
⋮----
atom = prev_config[idx]
⋮----
# 计算需要的z方向移动
surface_top_z = surface_relaxed.positions[:, 2].max()
ads_positions = np.array([pos for _, pos in ads_atoms_list])
ads_center_z = ads_positions[:, 2].mean()
z_shift = (surface_top_z + 3.5) - ads_center_z
⋮----
# 合并表面和移位后的分子
⋮----
new_pos = pos.copy()
⋮----
# 计算脱附态能量（不优化，因为这是gas phase）
⋮----
desorbed_result = self.predict_energy(
⋮----
relax=False,  # 不优化，保持gas phase几何
⋮----
# 脱附能 = E_desorbed - E_surface (应该是正值)
⋮----
# 保存脱附态结构
⋮----
# Fallback: use clean surface
⋮----
# 使用固定吸附位点
# 优先使用actual_reference_site（如果已从之前的吸附物获取）
# 这确保了所有吸附物使用完全相同的位点
site_to_use = actual_reference_site if actual_reference_site is not None else fixed_adsorption_site
⋮----
ads_result = self.predict_adsorption_energy(
⋮----
num_sites=1,  # 只试一个位点
⋮----
fixed_site_position=site_to_use,  # 传递固定位点
skip_surface_relaxation=True,  # 跳过表面优化（已经优化过了）
⋮----
# 正常模式 - 尝试多个吸附位点
⋮----
ads_name = str(ads)
⋮----
# 保存最佳配置
⋮----
# 更新参考吸附位点（从第一个成功处理的非"*"吸附物获取）
# 这个位点将用于所有后续吸附物，确保位点一致性
⋮----
actual_reference_site = binding_pos.copy()
⋮----
# 3. 计算相对能量（相对于最稳定的吸附质）
⋮----
detail = "; ".join([f"{name}: {err}" for name, err in failed_adsorbates]) or "unknown"
⋮----
min_energy = min(adsorbate_energies.values())
relative_energies = {
⋮----
def __del__(self)
⋮----
"""清理临时文件"""
⋮----
# 简单测试
````

## File: core/pathway/free_energy_router.py
````python
"""Free energy estimation routing for mechanism search.

Dispatches state-level energy queries to the appropriate backend
(thermal UMA, electrochemical CHE, heuristic) and returns unified
``StateEnergyEstimate``-compatible dicts.
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
# ---------------------------------------------------------------------------
# Backend implementations
⋮----
class ThermalStateEnergyBackend
⋮----
"""Estimate adsorption free energy via UMA calculator.

    Wraps the existing ``compute_adsorption_energies`` tool, but returns
    a simplified energy estimate.
    """
⋮----
def __init__(self, tools: Any = None)
⋮----
"""Return a dict matching StateEnergyEstimate fields."""
⋮----
call_kwargs: Dict[str, Any] = {
# Use fixed adsorption site if provided (from Agent2/3 best site)
fixed_site = kwargs.get("fixed_site")
⋮----
call_kwargs["num_sites"] = 1  # only one site when fixed
result = self._tools.compute_adsorption_energies(**call_kwargs)
e_ads = result.get("adsorption_energies", {}).get(species_label)
⋮----
@staticmethod
    def _heuristic_estimate(species_label: str) -> Dict[str, Any]
⋮----
"""Placeholder heuristic when UMA is unavailable."""
⋮----
class ElectroStateEnergyBackend
⋮----
"""Computational Hydrogen Electrode (CHE) free energy estimate.

    Phase 1 minimal implementation: applies CHE correction
    ``dG = dE + dZPE - TdS + neU + 0.0592*pH*n``
    using tabulated ZPE/entropy values where available.
    """
⋮----
"""Apply CHE correction to an adsorption energy."""
⋮----
# CHE correction: ΔG = ΔE + neU + 0.0592 * pH * n
che_correction = n_electrons * self.voltage_V + 0.0592 * self.pH * n_electrons
free_energy = adsorption_energy_eV + che_correction
⋮----
# Router
⋮----
class FreeEnergyRouter
⋮----
"""Routes free energy estimation to the appropriate backend."""
⋮----
"""Estimate free energy for a single state."""
backend = backend or self.default_backend
⋮----
"""Batch estimate free energies for multiple states."""
⋮----
__all__ = [
````

## File: core/pathway/improved_barrier_predictor.py
````python
"""
Improved NEB Barrier Predictor - 支持原子转移反应

这个模块解决了原始 barrier_predictor.py 的关键问题：
当反应物和产物原子数不同时，仍然能够计算能垒。

关键改进：
1. 智能识别需要添加/移除的原子
2. 自动在表面上放置额外的原子
3. 构建合理的NEB初始路径
4. 支持多种反应类型（吸附、解离、表面反应等）

作者: Claude
日期: 2026-01-21
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
@dataclass
class AtomTransferSetup
⋮----
"""原子转移反应的设置信息"""
reaction_type: str  # "add", "remove", "exchange"
atoms_to_add: List[str]  # 需要添加的原子
atoms_to_remove: List[int]  # 需要移除的原子索引
placement_strategy: str  # "nearby", "far", "specific"
placement_sites: List[np.ndarray]  # 放置位置
description: str  # 人类可读的描述
⋮----
class ImprovedBarrierPredictor
⋮----
"""
    改进的反应能垒预测器

    新功能：
    - 处理原子数不匹配的反应
    - 自动识别反应类型
    - 智能放置额外原子
    - 构建合理的NEB路径

    Example:
        CO + O → CO2 (原子转移反应):
        1. 识别需要一个O原子
        2. 在CO附近放置O
        3. 构建 [CO, O_nearby] → [CO2] 的NEB路径
    """
⋮----
"""
        分析反应物和产物的原子差异

        Parameters
        ----------
        reactant : Atoms
            反应物结构
        product : Atoms
            产物结构

        Returns
        -------
        AtomTransferSetup
            原子转移的设置信息
        """
r_symbols = reactant.get_chemical_symbols()
p_symbols = product.get_chemical_symbols()
⋮----
n_reactant = len(r_symbols)
n_product = len(p_symbols)
⋮----
# 产物比反应物多 → 需要添加原子
# 找出多了哪些原子
⋮----
r_count = Counter(r_symbols)
p_count = Counter(p_symbols)
⋮----
atoms_to_add = []
⋮----
diff = count - r_count.get(element, 0)
⋮----
placement_sites=[],  # 将在后续计算
⋮----
# 产物比反应物少 → 需要移除原子
⋮----
atoms_to_remove_symbols = []
⋮----
diff = count - p_count.get(element, 0)
⋮----
atoms_to_remove=[],  # 将在后续计算具体索引
⋮----
def find_surface_atoms(self, structure: Atoms, z_threshold: float = 0.5) -> List[int]
⋮----
"""
        找出表面原子（最上层）

        Parameters
        ----------
        structure : Atoms
            表面+吸附质结构
        z_threshold : float
            判断表面的z坐标阈值（相对于最高点，单位Å）

        Returns
        -------
        List[int]
            表面原子的索引列表
        """
positions = structure.get_positions()
z_coords = positions[:, 2]
z_max = z_coords.max()
⋮----
surface_indices = [
⋮----
"""
        找出吸附质原子（不是表面的原子）

        Parameters
        ----------
        structure : Atoms
            表面+吸附质结构
        min_z_above_surface : float
            吸附质距离表面的最小高度（Å）

        Returns
        -------
        List[int]
            吸附质原子的索引列表
        """
⋮----
# 找到表面的平均高度
# 假设最下面50%的原子是表面/slab
z_sorted = sorted(z_coords)
median_idx = len(z_sorted) // 2
z_surface = z_sorted[median_idx]
⋮----
adsorbate_indices = [
⋮----
"""
        在吸附质附近的表面上放置一个原子

        Parameters
        ----------
        structure : Atoms
            原始结构
        adsorbate_indices : List[int]
            吸附质原子的索引
        atom_symbol : str
            要放置的原子符号（如 "O", "H"）
        distance : float
            与吸附质的距离（Å）
        height_above_surface : float
            距离表面的高度（Å）

        Returns
        -------
        new_structure : Atoms
            添加了原子后的结构
        new_atom_index : int
            新添加原子的索引
        """
# 计算吸附质的中心位置
ads_positions = structure.get_positions()[adsorbate_indices]
ads_center = ads_positions.mean(axis=0)
⋮----
# 找到表面高度
all_positions = structure.get_positions()
z_coords = all_positions[:, 2]
surface_z = np.median(z_coords)
⋮----
# 在吸附质旁边放置新原子
# 方向：x方向偏移
new_pos = ads_center.copy()
new_pos[0] += distance  # x方向偏移
new_pos[2] = surface_z + height_above_surface  # z高度
⋮----
# 添加原子
new_structure = structure.copy()
⋮----
new_atom_index = len(new_structure) - 1
⋮----
"""
        为涉及原子转移的反应构建NEB路径

        Strategy:
        1. 分析原子差异
        2. 如果需要添加原子，在反应物中添加
        3. 构建从 (reactant + extra atoms) 到 product 的路径

        Parameters
        ----------
        reactant : Atoms
            反应物结构
        product : Atoms
            产物结构
        n_frames : int
            NEB帧数

        Returns
        -------
        List[Atoms]
            NEB路径（所有帧原子数相同）
        """
# 分析原子差异
transfer_setup = self.analyze_atom_difference(reactant, product)
⋮----
# 原子数相同，直接插值
⋮----
# 需要添加原子
⋮----
# 找到吸附质原子
ads_indices = self.find_adsorbate_atoms(reactant)
⋮----
# 添加每个需要的原子
modified_reactant = reactant.copy()
⋮----
distance=3.0,  # 3Å 距离
⋮----
# 现在原子数应该匹配
⋮----
# 插值
⋮----
# 需要移除原子（脱附）
⋮----
# 在产物中添加远离的原子（模拟脱附）
# 这里简化处理：将脱附的原子移到远处
modified_product = product.copy()
⋮----
# 找出反应物中哪些原子在产物中不存在
# （这里简化：假设是最后几个原子）
n_to_remove = len(reactant) - len(product)
⋮----
# 在产物上方远处添加原子
pos = modified_product.get_positions().mean(axis=0)
pos[2] += 10.0  # 10Å 上方
⋮----
# 找出需要移除的原子类型
removed_symbol = reactant.get_chemical_symbols()[-1]
⋮----
"""
        简单的线性插值

        Parameters
        ----------
        initial : Atoms
            初始结构
        final : Atoms
            最终结构
        n_frames : int
            帧数

        Returns
        -------
        List[Atoms]
            插值后的路径
        """
⋮----
frames = [initial]
⋮----
frame = initial.copy()
⋮----
# ASE 插值
⋮----
"""
        从反应物和产物结构预测能垒（支持原子数不同）

        Parameters
        ----------
        reactant : str or Atoms
            反应物结构
        product : str or Atoms
            产物结构
        n_frames : int
            NEB 帧数
        fmax : float
            力收敛标准
        max_steps : int
            最大优化步数
        reaction_name : str, optional
            反应名称

        Returns
        -------
        Dict
            包含能垒和反应能的结果
        """
# 读取结构
⋮----
reactant_atoms = read(reactant)
⋮----
reactant_atoms = reactant.copy()
⋮----
product_atoms = read(product)
⋮----
product_atoms = product.copy()
⋮----
# 构建NEB路径（自动处理原子数不匹配）
⋮----
frames = self.build_neb_path_with_atom_transfer(
⋮----
# 运行 NEB (这里简化，实际应该调用 Fairchem)
# TODO: 实现实际的NEB优化
⋮----
def test_atom_transfer_neb()
⋮----
"""测试原子转移NEB"""
⋮----
# 创建一个简单的Pt(111)表面
slab = fcc111('Pt', size=(3, 3, 4), vacuum=10.0)
⋮----
# 反应物: CO adsorbed
co = molecule('CO')
⋮----
reactant = slab.copy()
⋮----
# 产物: CO2 adsorbed (比反应物多1个O原子)
co2 = molecule('CO2')
slab_for_product = fcc111('Pt', size=(3, 3, 4), vacuum=10.0)
⋮----
product = slab_for_product.copy()
⋮----
# 测试改进的预测器
predictor = ImprovedBarrierPredictor(
⋮----
result = predictor.predict_from_structures(
````

## File: core/pathway/initial_final_validator.py
````python
"""
初末态结构验证器 (Agent 4.5)

自动验证和修复NEB计算的初始态和最终态结构，确保：
1. 键长合理
2. 无原子重叠
3. 配位数正确
4. 吸附位点合理
5. 氢供体位置合适

支持任意反应、任意体系、任意吸附分子的通用性验证。
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
class ReactionRuleDatabase
⋮----
"""
    反应类型规则数据库
    """
⋮----
# 标准共价键长（Å）
STANDARD_BOND_LENGTHS = {
⋮----
'C-surface': 2.0,  # 典型C吸附在金属表面
'H-surface': 1.8,  # 典型H吸附在金属表面
'O-surface': 2.0,  # 典型O吸附在金属表面
⋮----
# 允许的键长范围（倍数）
BOND_LENGTH_TOLERANCE = 0.3  # ±30%
⋮----
# 反应类型模板
REACTION_TEMPLATES = {
⋮----
"C-H_forming": (1.0, 2.5),  # 形成中的C-H键
"H-surface": (1.5, 2.2),    # H供体到表面距离
"C-surface": (1.8, 2.5),    # C到表面距离
⋮----
"C": {"initial": (3, 4), "final": (4, 4)},  # sp2 → sp3
⋮----
"C-H_breaking": (1.0, 1.3),  # 初态的C-H键
⋮----
"C": {"initial": (4, 4), "final": (3, 4)},  # sp3 → sp2
⋮----
"C-O_forming": (1.2, 2.8),  # 形成中的C-O键
⋮----
"min_interatomic_distance": 0.8,  # 任意两原子最小距离
"max_bond_length_factor": 1.5,   # 共价键最大长度倍数
⋮----
@classmethod
    def get_rule(cls, reaction_type: str) -> Dict
⋮----
"""获取反应类型的规则"""
⋮----
@classmethod
    def infer_reaction_type(cls, initial_formula: str, final_formula: str) -> str
⋮----
"""
        从化学式推断反应类型
        
        Examples:
            "*CH2", "*CH3" → "hydrogenation"
            "*CO + *O", "*CO2" → "c_o_coupling"
        """
# 简单的启发式推断
initial_h = initial_formula.count('H')
final_h = final_formula.count('H')
⋮----
class StructureValidator
⋮----
"""
    结构合理性验证器
    """
⋮----
def __init__(self, reaction_type: str = "generic")
⋮----
"""
        检查键长是否在合理范围
        
        Args:
            atoms: ASE Atoms对象
            bond_pairs: 键的原子索引对列表 [(i, j), ...]
            bond_type: 键类型（如 "C-H_forming"）
        
        Returns:
            (is_valid, message)
        """
⋮----
distance = atoms.get_distance(i, j)
⋮----
def check_atom_overlap(self, atoms: Atoms, min_distance: float = 0.8) -> Tuple[bool, str]
⋮----
"""
        检查原子是否过于接近（重叠）
        
        Args:
            atoms: ASE Atoms对象
            min_distance: 最小允许距离（Å）
        
        Returns:
            (is_valid, message)
        """
positions = atoms.get_positions()
⋮----
# 计算所有原子对的距离
distances = cdist(positions, positions)
np.fill_diagonal(distances, np.inf)  # 忽略自身
⋮----
min_dist = distances.min()
overlap_indices = np.where(distances == min_dist)
⋮----
"""
        检查原子的配位数
        
        Args:
            atoms: ASE Atoms对象
            atom_indices: 要检查的原子索引列表
            expected_range: 期望的配位数范围 (min, max)
            cutoff: 配位数统计的距离阈值（Å）
        
        Returns:
            (is_valid, message)
        """
⋮----
element = atoms[idx].symbol
⋮----
# 计算近邻数
⋮----
distances = np.linalg.norm(positions - positions[idx], axis=1)
neighbors = np.sum((distances > 0.1) & (distances < cutoff))
⋮----
"""
        检查氢供体到碳受体的距离（加氢反应专用）
        
        Args:
            atoms: ASE Atoms对象
            h_donor_idx: 氢供体原子索引
            c_acceptor_idx: 碳受体原子索引
        
        Returns:
            (is_valid, message)
        """
⋮----
distance = atoms.get_distance(h_donor_idx, c_acceptor_idx)
⋮----
"""
        检查吸附物到表面的距离
        
        Args:
            atoms: ASE Atoms对象
            adsorbate_indices: 吸附物原子索引
            surface_indices: 表面原子索引
            min_distance: 最小距离（太近可能重叠）
            max_distance: 最大距离（太远无法吸附）
        
        Returns:
            (is_valid, message)
        """
adsorbate_pos = atoms.get_positions()[adsorbate_indices]
surface_pos = atoms.get_positions()[surface_indices]
⋮----
# 计算吸附物到表面的最小距离
distances = cdist(adsorbate_pos, surface_pos)
⋮----
class StructureFixer
⋮----
"""
    结构自动修复工具
    """
⋮----
def __init__(self)
⋮----
"""
        调整两个原子间的距离到目标值
        
        Args:
            atoms: ASE Atoms对象
            atom_i, atom_j: 两个原子的索引
            target_distance: 目标距离（Å）
            freeze_indices: 冻结的原子索引（通常是表面原子）
        
        Returns:
            修复后的Atoms对象
        """
atoms = atoms.copy()
⋮----
# 确定哪个原子可移动
⋮----
# 都可移动，移动第二个
⋮----
# 计算方向向量
direction = atoms[movable_idx].position - atoms[anchor_idx].position
current_distance = np.linalg.norm(direction)
⋮----
# 设置新位置
new_position = atoms[anchor_idx].position + direction * target_distance
⋮----
"""
        分离重叠的原子
        
        Args:
            atoms: ASE Atoms对象
            atom_i, atom_j: 重叠的两个原子
            min_distance: 最小安全距离
            freeze_indices: 冻结的原子索引
        
        Returns:
            修复后的Atoms对象
        """
⋮----
"""
        调整吸附物到表面的高度（Z方向）
        
        Args:
            atoms: ASE Atoms对象
            adsorbate_indices: 吸附物原子索引
            target_height: 目标高度（Å，相对于表面顶层）
        
        Returns:
            修复后的Atoms对象
        """
⋮----
# 找到表面顶层Z坐标
all_z = atoms.get_positions()[:, 2]
adsorbate_z = all_z[adsorbate_indices]
surface_z = np.delete(all_z, adsorbate_indices)
surface_top = surface_z.max()
⋮----
# 计算吸附物质心当前高度
current_height = adsorbate_z.mean() - surface_top
⋮----
# 整体移动吸附物
shift = target_height - current_height
⋮----
class InitialFinalStateValidator
⋮----
"""
    Agent 4.5: 初末态结构验证与修复
    
    工作流程:
    1. 接收反应初末态结构
    2. 运行一系列验证检查
    3. 对于失败的检查，尝试自动修复
    4. 生成验证报告
    5. 返回修复后的结构或标记为需要人工审查
    """
⋮----
def __init__(self, reaction_type: str = "generic", auto_fix: bool = True)
⋮----
"""
        Args:
            reaction_type: 反应类型（hydrogenation, dehydrogenation, c_o_coupling, generic）
            auto_fix: 是否自动修复问题
        """
⋮----
"""
        验证并修复初末态结构
        
        Args:
            initial_state: 初始态结构
            final_state: 最终态结构
            adsorbate_indices: 吸附物原子索引
            surface_indices: 表面原子索引（冻结）
            reaction_info: 反应的额外信息（如氢供体位置等）
        
        Returns:
            (fixed_initial, fixed_final, validation_report)
        """
⋮----
# 初始化报告
report = {
⋮----
# 验证初始态
⋮----
# 验证最终态
⋮----
# 判断整体状态
has_errors = (len(initial_report["errors"]) > 0) or (len(final_report["errors"]) > 0)
⋮----
# 打印摘要
⋮----
"""
        验证单个状态（初态或末态）
        
        Returns:
            (fixed_atoms, report_dict)
        """
⋮----
# 检查1: 原子重叠
⋮----
# 提取重叠的原子索引
# 这里简化处理，实际需要解析msg
⋮----
# TODO: 实现智能分离
⋮----
# 检查2: 表面距离
⋮----
atoms = self.fixer.adjust_adsorbate_height(atoms, adsorbate_indices)
# 重新检查
⋮----
# 检查3: 氢供体距离（如果是加氢反应）
⋮----
target_distance = 2.0  # 合理的H-C距离
atoms = self.fixer.fix_bond_distance(
⋮----
# 检查4: 配位数（根据反应类型）
⋮----
coord_rules = self.validator.rules["coordination_changes"].get("C", None)
⋮----
# 找到碳原子索引（简化：假设第一个吸附原子是C）
c_indices = [i for i in adsorbate_indices if atoms[i].symbol == 'C']
⋮----
expected = coord_rules.get("initial", (3, 4))
⋮----
expected = coord_rules.get("final", (3, 4))
⋮----
def _print_report(self, report: Dict)
⋮----
"""打印验证报告"""
⋮----
state_name = state.replace("_state", "").upper()
⋮----
for p in report[state]["passed"][:3]:  # 只显示前3个
⋮----
def quick_test()
⋮----
"""快速测试"""
⋮----
# 创建测试表面
slab = fcc111('Pt', size=(3, 3, 4), vacuum=10.0)
⋮----
# 添加*CH2
⋮----
c_idx = len(slab) - 1
⋮----
# 添加2个H（模拟CH2）
c_pos = slab[c_idx].position
⋮----
initial_state = slab.copy()
⋮----
# 创建最终态（添加第3个H，模拟*CH3）
# 故意让H原子距离很远，测试验证器
final_state = slab.copy()
final_state.append(Atom('H', position=c_pos + [0, 0, 4.5]))  # 太远！
⋮----
# 定义吸附物和表面索引
n_surface = len(fcc111('Pt', size=(3, 3, 4), vacuum=10.0))
adsorbate_indices = list(range(n_surface, len(initial_state)))
surface_indices = list(range(n_surface))
⋮----
# 运行验证
validator = InitialFinalStateValidator(reaction_type="hydrogenation", auto_fix=True)
⋮----
"h_donor_idx": len(final_state) - 1,  # 最后一个H
````

## File: core/pathway/llm_neb_controller.py
````python
"""
LLM-Controlled NEB Structure Preparation

This module provides a general-purpose, LLM-controlled approach for preparing
NEB (Nudged Elastic Band) calculations for any reaction on any surface.

The LLM fully controls:
1. Analysis of reactant and product structures
2. Decision on whether atoms need to be added/removed
3. Positions of added atoms
4. Validation of prepared structures
5. Iterative refinement if structures are unreasonable

Key principle: The adsorption sites are pre-determined and cannot be changed.
LLM only decides how to match atoms between reactant and product for NEB.
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
@dataclass
class NEBStructurePlan
⋮----
"""LLM's plan for NEB structure preparation"""
reaction_type: str  # e.g., "hydrogenation", "dehydrogenation", "bond_breaking"
atoms_to_add_to_reactant: List[Dict]  # [{"element": "H", "position": [x,y,z], "reason": "..."}]
atoms_to_add_to_product: List[Dict]
atoms_to_remove_from_reactant: List[int]  # atom indices
atoms_to_remove_from_product: List[int]
atom_mapping: Dict[int, int]  # reactant_idx -> product_idx for alignment
reasoning: str
confidence: float  # 0-1
⋮----
# Energy correction information
original_reactant_energy: Optional[float] = None
adjusted_reactant_energy: Optional[float] = None
original_product_energy: Optional[float] = None
adjusted_product_energy: Optional[float] = None
⋮----
@dataclass
class ValidationResult
⋮----
"""LLM's validation of prepared structures"""
is_valid: bool
issues: List[str]
suggestions: List[str]
corrected_plan: Optional[NEBStructurePlan] = None
auto_fixes_applied: List[str] = None  # Record what was auto-fixed
⋮----
def __post_init__(self)
⋮----
class LLMNEBController
⋮----
"""
    LLM-controlled NEB structure preparation.

    This class uses an LLM to intelligently prepare structures for NEB calculations
    by analyzing the chemistry and deciding how to adjust atom counts and positions.
    """
⋮----
"""
        Initialize the LLM NEB Controller.

        Parameters
        ----------
        model : str, optional
            LLM model to use. Defaults to environment variable or claude-sonnet.
        max_retries : int
            Maximum attempts for structure preparation
        temperature : float
            LLM temperature for responses
        """
⋮----
def _init_llm(self)
⋮----
"""Initialize LLM client"""
⋮----
api_key = os.getenv("OPENAI_API_KEY")
base_url = os.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1")
⋮----
"""
        Convert ASE Atoms to a detailed string representation for LLM.

        Includes:
        - Cell parameters
        - All atom positions with indices and markers
        - Adsorbate atoms clearly marked
        - Interatomic distances for adsorbate atoms
        """
lines = [f"=== {name} ==="]
⋮----
# Identify surface element(s)
symbols = atoms.get_chemical_symbols()
⋮----
symbol_counts = Counter(symbols)
surface_element = max(symbol_counts, key=symbol_counts.get)
⋮----
# Cell
cell = atoms.cell
⋮----
# Surface z-range
positions = atoms.get_positions()
surface_z = [positions[i, 2] for i in range(len(atoms)) if symbols[i] == surface_element]
⋮----
# Adsorbate atoms (detailed)
⋮----
pos = positions[idx]
sym = symbols[idx]
⋮----
# Distances between adsorbate atoms
⋮----
dist = np.linalg.norm(positions[idx1] - positions[idx2])
⋮----
# Full coordinates (POSCAR-like)
⋮----
note = "ADSORBATE" if i in adsorbate_indices else ""
⋮----
def _call_llm(self, system_prompt: str, user_prompt: str) -> str
⋮----
"""Call LLM and return response"""
⋮----
response = self._client.chat.completions.create(
⋮----
"""
        Use LLM to analyze the reaction and create a plan for NEB structure preparation.

        Parameters
        ----------
        reactant : Atoms
            Reactant structure (already relaxed at its adsorption site)
        product : Atoms
            Product structure (already relaxed at its adsorption site)
        reactant_adsorbate_indices : List[int]
            Indices of adsorbate atoms in reactant
        product_adsorbate_indices : List[int]
            Indices of adsorbate atoms in product
        reaction_name : str, optional
            Human-readable reaction name

        Returns
        -------
        NEBStructurePlan
            LLM's plan for structure preparation
        """
# Prepare detailed structure descriptions
reactant_desc = self._atoms_to_detailed_string(
product_desc = self._atoms_to_detailed_string(
⋮----
# Analyze atom count differences
r_symbols = reactant.get_chemical_symbols()
p_symbols = product.get_chemical_symbols()
⋮----
r_counts = Counter(r_symbols)
p_counts = Counter(p_symbols)
⋮----
all_elements = set(r_counts.keys()) | set(p_counts.keys())
diff_info = []
⋮----
r_n = r_counts.get(elem, 0)
p_n = p_counts.get(elem, 0)
⋮----
diff_str = "\n".join(diff_info) if diff_info else "  No difference in atom counts"
⋮----
system_prompt = """You are an expert computational chemist specializing in catalytic reactions and NEB (Nudged Elastic Band) calculations.
⋮----
# Get current adsorbate sequences
r_ads_symbols = [reactant[i].symbol for i in reactant_adsorbate_indices]
p_ads_symbols = [product[i].symbol for i in product_adsorbate_indices]
⋮----
# Prepare pathway continuity information if provided
pathway_context = ""
⋮----
prev_ads_symbols = [previous_step_product[i].symbol for i in previous_step_product_adsorbate_indices]
prev_ads_positions = previous_step_product.positions[previous_step_product_adsorbate_indices]
⋮----
# Describe previous step's product
pathway_context = f"""
⋮----
user_prompt = f"""Analyze this surface catalytic reaction and create a NEB structure preparation plan.
⋮----
response = self._call_llm(system_prompt, user_prompt)
⋮----
# Parse JSON from response
json_match = re.search(r'```json\s*(.*?)\s*```', response, re.DOTALL)
⋮----
json_match = re.search(r'\{.*\}', response, re.DOTALL)
⋮----
plan_dict = json.loads(json_match.group(1) if json_match.group(1) else json_match.group())
⋮----
plan_dict = json.loads(json_match.group())
⋮----
# Convert to NEBStructurePlan
plan = NEBStructurePlan(
⋮----
"""
        Apply the LLM's plan to prepare structures for NEB.

        Returns
        -------
        Tuple of (adjusted_reactant, adjusted_product, new_reactant_indices, new_product_indices)
        """
# Make copies
adj_reactant = reactant.copy()
adj_product = product.copy()
new_r_indices = list(reactant_adsorbate_indices)
new_p_indices = list(product_adsorbate_indices)
⋮----
# Remove atoms first (if any)
⋮----
indices_to_remove = sorted(plan.atoms_to_remove_from_reactant, reverse=True)
⋮----
new_r_indices = [i if i < idx else i - 1 for i in new_r_indices if i != idx]
⋮----
indices_to_remove = sorted(plan.atoms_to_remove_from_product, reverse=True)
⋮----
new_p_indices = [i if i < idx else i - 1 for i in new_p_indices if i != idx]
⋮----
# Add atoms
from ase import Atom  # Import once at the top of the section
⋮----
element = atom_info["element"]
position = np.array(atom_info["position"])
insert_index = atom_info.get("insert_index", None)
⋮----
# Validate position
positions = adj_reactant.get_positions()
⋮----
min_dist = np.min(np.linalg.norm(positions - position, axis=1))
⋮----
# Find nearest atom and move away
nearest_idx = np.argmin(np.linalg.norm(positions - position, axis=1))
direction = position - positions[nearest_idx]
direction = direction / np.linalg.norm(direction)
position = positions[nearest_idx] + direction * 1.0
⋮----
# Insert atom at specified index or append
⋮----
# Insert at specific position
⋮----
new_atom = Atom(element, position=position)
temp_atoms = adj_reactant[:insert_index] + AseAtoms([new_atom]) + adj_reactant[insert_index:]
adj_reactant = temp_atoms
⋮----
# Update adsorbate indices: indices >= insert_index shift by 1
new_r_indices = [i if i < insert_index else i + 1 for i in new_r_indices]
⋮----
# Insert new atom index at correct position (maintain sorted order)
insert_pos = sum(1 for idx in new_r_indices if idx < insert_index)
⋮----
# Append at end (default behavior)
⋮----
positions = adj_product.get_positions()
⋮----
temp_atoms = adj_product[:insert_index] + AseAtoms([new_atom]) + adj_product[insert_index:]
adj_product = temp_atoms
⋮----
new_p_indices = [i if i < insert_index else i + 1 for i in new_p_indices]
⋮----
insert_pos = sum(1 for idx in new_p_indices if idx < insert_index)
⋮----
# Step: Handle atom reordering to match product atom types
⋮----
# Check if atom ordering matches by element type
reactant_symbols = [adj_reactant[i].symbol for i in new_r_indices]
product_symbols = [adj_product[i].symbol for i in new_p_indices]
⋮----
# Simple reordering strategy: match atoms by element type
# For each position in product, find a reactant atom of the same element
new_r_indices_reordered = []
used_reactant_indices = set()
⋮----
product_elem = adj_product[p_idx].symbol
⋮----
# Find an unused reactant atom with the same element
found = False
⋮----
found = True
⋮----
# ⭐ CRITICAL FIX: Actually reorder the reactant atoms physically
⋮----
# Step 1: Identify which atoms are surface vs adsorbate
n_atoms = len(adj_reactant)
n_ads = len(new_r_indices)
adsorbate_indices_set = set(new_r_indices)
⋮----
# Step 2: Build complete reordering mapping
# Surface atoms: keep in original order (indices 0 to n-n_ads-1)
# Adsorbate atoms: reorder according to new_r_indices_reordered
reorder_sequence = []
⋮----
# First, add all surface atoms (not in adsorbate_indices)
⋮----
# Then, add reordered adsorbate atoms
⋮----
# Step 3: Apply reordering using ASE slicing
adj_reactant = adj_reactant[reorder_sequence]
⋮----
# Step 4: Update adsorbate indices (they're now at the end)
n_surface_atoms = len(reorder_sequence) - n_ads
new_r_indices = list(range(n_surface_atoms, len(adj_reactant)))
⋮----
# Verify the reordering
reactant_symbols_new = [adj_reactant[i].symbol for i in new_r_indices]
⋮----
"""
        Automatically fix structural issues based on LLM validation results.

        This is the key agent loop: LLM detects problems → LLM proposes fixes → apply fixes

        The LLM sees:
        - The problematic structures
        - Validation issues identified
        - Validation suggestions from previous LLM

        And returns:
        - Corrected atom positions
        - Reordering instructions

        Returns
        -------
        Tuple of (fixed_reactant, fixed_product, new_r_indices, new_p_indices, fixes_applied)
        """
⋮----
# Prepare detailed description of the problem
⋮----
issues_str = "\n".join(f"  - {issue}" for issue in validation.issues)
suggestions_str = "\n".join(f"  - {sug}" for sug in validation.suggestions)
⋮----
system_prompt = """You are an expert computational chemist fixing NEB structure problems.
⋮----
user_prompt = f"""Fix the following NEB structure problems:
⋮----
# Parse JSON
⋮----
json_match = re.search(r'\{.*"fixes".*\}', response, re.DOTALL)
⋮----
fix_dict = json.loads(json_match.group(1) if '```' in response else json_match.group())
⋮----
# Apply LLM-generated fixes
⋮----
"""Apply fixes generated by LLM"""
fixes_applied = []
⋮----
new_r_idx = list(reactant_adsorbate_indices)
new_p_idx = list(product_adsorbate_indices)
⋮----
fixes = fix_dict.get("fixes", [])
⋮----
fix_type = fix.get("type")
description = fix.get("description", "")
⋮----
# Reorder atoms according to mapping
mapping = fix.get("reactant_atom_mapping", {})
⋮----
# Convert string keys to int
swap_pairs = []
⋮----
old_idx = int(old_str)
new_idx = int(new_str) if isinstance(new_str, (int, str)) else new_str
⋮----
new_idx = int(new_idx)
⋮----
# Create a new Atoms object with swapped atoms
# This ensures the swap is actually applied
⋮----
# Save original tags before rebuilding
original_tags = adj_reactant.get_tags()
⋮----
new_atoms_list = []
⋮----
# Check if this index should be swapped
target_idx = i
⋮----
target_idx = new_idx
⋮----
target_idx = old_idx
⋮----
# Copy the atom from target index
atom = adj_reactant[target_idx]
⋮----
# Rebuild structure with TAGS preserved
adj_reactant = AseAtoms(new_atoms_list,
⋮----
# Apply swapped tags (tags should follow the atoms)
new_tags = []
⋮----
structure = fix.get("structure", "reactant")
atom_idx = int(fix.get("atom_index"))
new_pos = np.array(fix.get("new_position"))
reason = fix.get("reason", "")
⋮----
"""
        Rule-based fallback fixes (original implementation)
        """
⋮----
# Fix 1: Atom ordering mismatch
# Check if adsorbate atoms have different elements at same indices
r_symbols = np.array(adj_reactant.get_chemical_symbols())
p_symbols = np.array(adj_product.get_chemical_symbols())
⋮----
# Find mismatches in adsorbate region
mismatches = []
⋮----
# Strategy: Reorder reactant adsorbate atoms to match product
# Build a mapping based on element types
r_ads_positions = adj_reactant.get_positions()[new_r_idx]
p_ads_positions = adj_product.get_positions()[new_p_idx]
r_ads_symbols = [r_symbols[i] for i in new_r_idx]
p_ads_symbols = [p_symbols[i] for i in new_p_idx]
⋮----
# For each product adsorbate atom, find the closest same-element reactant atom
⋮----
n_ads = len(new_r_idx)
cost_matrix = np.full((n_ads, n_ads), 1e10)
⋮----
# Cost is distance between atoms
⋮----
# Reorder reactant adsorbate atoms
new_order_r_idx = [new_r_idx[col_ind[i]] for i in range(n_ads)]
⋮----
# Create a full reordering array for the entire reactant structure
# Keep surface atoms in place, only reorder adsorbate
surface_indices = [i for i in range(len(adj_reactant)) if i not in new_r_idx]
full_new_order = surface_indices + new_order_r_idx
⋮----
adj_reactant = adj_reactant[full_new_order]
⋮----
# Update indices
new_r_idx = list(range(len(surface_indices), len(adj_reactant)))
⋮----
# Fix 2: Atoms too close (< 0.8 Å)
r_pos = adj_reactant.get_positions()
p_pos = adj_product.get_positions()
⋮----
# Check reactant adsorbate atoms
⋮----
dist = np.linalg.norm(r_pos[idx_i] - r_pos[idx_j])
⋮----
# Move atoms apart along the line connecting them
direction = r_pos[idx_i] - r_pos[idx_j]
⋮----
# Move each atom 0.5 Å away from midpoint
midpoint = (r_pos[idx_i] + r_pos[idx_j]) / 2
⋮----
# Check product adsorbate atoms
⋮----
dist = np.linalg.norm(p_pos[idx_i] - p_pos[idx_j])
⋮----
direction = p_pos[idx_i] - p_pos[idx_j]
⋮----
midpoint = (p_pos[idx_i] + p_pos[idx_j]) / 2
⋮----
# Fix 3: For hydrogenation, ensure H is approaching from reasonable distance
# If reaction type is hydrogenation and H is too close to target
⋮----
# We added atoms, likely hydrogenation
added_indices = [i for i in new_r_idx if i >= len(reactant)]
⋮----
# Find the closest non-H atom in adsorbate
other_ads_indices = [i for i in new_r_idx if i != added_idx and r_symbols[i] != 'H']
⋮----
distances = [np.linalg.norm(r_pos[added_idx] - r_pos[other_idx])
min_dist = min(distances)
closest_idx = other_ads_indices[np.argmin(distances)]
⋮----
# Move H further away (to ~2.0-2.5 Å)
target_dist = 2.2
direction = r_pos[added_idx] - r_pos[closest_idx]
⋮----
new_h_pos = r_pos[closest_idx] + direction * target_dist
⋮----
"""
        Use LLM to validate the prepared structures.

        Returns
        -------
        ValidationResult
            LLM's assessment of structure validity
        """
# Basic checks first
issues = []
⋮----
# Check atom counts
⋮----
# Check for overlapping atoms
r_pos = reactant.get_positions()
p_pos = product.get_positions()
⋮----
dist = np.linalg.norm(r_pos[i] - r_pos[j])
⋮----
dist = np.linalg.norm(p_pos[i] - p_pos[j])
⋮----
# If basic checks fail, no need for LLM
⋮----
# LLM validation
⋮----
system_prompt = """You are an expert computational chemist validating structures for NEB calculation.
⋮----
user_prompt = f"""Validate these structures for NEB calculation.
⋮----
result_dict = json.loads(json_match.group(1) if '```' in json_match.group() else json_match.group())
⋮----
# Fallback: assume valid if no JSON parsing issues and basic checks passed
⋮----
"""
        Main method: prepare structures for NEB calculation.

        This method orchestrates the full LLM-controlled agent workflow:
        1. Analyze reaction and create plan
        2. Apply plan to structures
        3. Validate results
        4. AUTO-FIX issues if validation fails
        5. Re-validate and retry if needed

        Parameters
        ----------
        reactant, product : Atoms
            Input structures (pre-relaxed at their adsorption sites)
        reactant_adsorbate_indices, product_adsorbate_indices : List[int]
            Adsorbate atom indices
        reaction_name : str, optional
            Human-readable reaction name
        previous_step_product : Atoms, optional
            Product structure from previous step in pathway (for continuity)
        previous_step_product_adsorbate_indices : List[int], optional
            Adsorbate indices for previous step's product
        step_index : int, optional
            Index of this step in the pathway (for context)

        Returns
        -------
        Tuple of (adjusted_reactant, adjusted_product, reactant_indices, product_indices, plan)
        """
⋮----
adj_reactant = None
adj_product = None
new_r_idx = None
new_p_idx = None
plan = None
all_fixes = []
⋮----
# Step 1: Analyze and plan
plan = self.analyze_reaction_and_plan(
⋮----
# Step 2: Apply plan
⋮----
# Record original and adjusted structures for energy correction
# This allows us to compute the energy offset from virtual atoms later
plan.original_reactant_energy = len(reactant)  # Store atom count as proxy
⋮----
# Step 3: Validate
validation = self.validate_structures(
⋮----
# Step 4: AUTO-FIX - This is the key agent self-healing loop!
⋮----
# Step 5: Re-validate after fixes
validation2 = self.validate_structures(
⋮----
# Continue to next attempt with fresh plan
⋮----
# If all retries failed, do final auto-fix and return
⋮----
# One more auto-fix pass with more aggressive settings
⋮----
# Ensure atom counts match
⋮----
"""
        Align product atoms to match reactant ordering for NEB.

        This ensures smooth interpolation by matching corresponding atoms.
        """
⋮----
r_sym = reactant.get_chemical_symbols()
p_sym = product.get_chemical_symbols()
⋮----
# Build cost matrix for Hungarian algorithm
n = len(reactant)
cost_matrix = np.full((n, n), 1e10)
⋮----
# Solve assignment
⋮----
# Reorder product
aligned_product = product[col_ind]
⋮----
# Update adsorbate indices
new_p_indices = []
⋮----
new_idx = np.where(col_ind == old_idx)[0]
````

## File: core/pathway/neb_frame_prep.py
````python
"""Shared endpoint preparation for interpolation-sensitive NEB workflows."""
⋮----
@dataclass
class PreparedEndpointPair
⋮----
reactant: Atoms
product: Atoms
reactant_adsorbate_indices: List[int]
product_adsorbate_indices: List[int]
reordered: bool
⋮----
"""Reorder product atoms to match reactant atom order when symbols match."""
product_adsorbate_indices = list(product_adsorbate_indices or [])
⋮----
reactant_symbols = reactant.get_chemical_symbols()
product_symbols = product.get_chemical_symbols()
⋮----
copied = product.copy()
⋮----
reactant_positions = reactant.get_positions()
product_positions = product.get_positions()
⋮----
from scipy.optimize import linear_sum_assignment  # type: ignore
⋮----
linear_sum_assignment = None
⋮----
assignment: Dict[int, int] = {}
grouped: Dict[str, Tuple[List[int], List[int]]] = {}
⋮----
cost = np.zeros((len(reactant_group), len(product_group)), dtype=float)
⋮----
remaining = list(product_group)
⋮----
chosen = min(
⋮----
reorder_sequence = [assignment[idx] for idx in range(len(reactant_symbols))]
⋮----
reordered_product = product[reorder_sequence]
old_to_new = {old: new for new, old in enumerate(reorder_sequence)}
new_ads = sorted(old_to_new[idx] for idx in product_adsorbate_indices if idx in old_to_new)
new_ads_set = set(new_ads)
⋮----
"""Prepare endpoints so interpolation follows the nearest periodic images per fragment."""
prepared_reactant = reactant.copy()
prepared_product = product.copy()
⋮----
reactant_ads = _valid_indices(
product_ads = _valid_indices(
⋮----
reordered = False
⋮----
reactant_graph = _build_adsorbate_graph(prepared_reactant, reactant_ads)
product_graph = _build_adsorbate_graph(prepared_product, product_ads)
⋮----
reactant_components = _connected_components(reactant_ads, reactant_graph)
product_components = _connected_components(product_ads, product_graph)
⋮----
ads_set_product = set(product_ads)
⋮----
nearest = _nearest_image_position(
⋮----
product_symbols = prepared_product.get_chemical_symbols()
⋮----
anchor_idx = _choose_component_anchor(component, product_symbols)
nearest_anchor = _nearest_image_position(
translation = nearest_anchor - prepared_product.positions[anchor_idx]
⋮----
"""Compute simple linear-interpolation distance diagnostics on prepared endpoints."""
⋮----
adsorbate_indices = _valid_indices(reactant, adsorbate_indices or [])
frame_count = max(3, int(n_images))
reactant_pos = np.array(reactant.get_positions(), dtype=float)
product_pos = np.array(product.get_positions(), dtype=float)
frame_atoms = reactant.copy()
⋮----
min_any_pair = float("inf")
min_any_pair_frame = -1
min_ads_surface = float("inf")
min_ads_surface_frame = -1
⋮----
alpha = frame_idx / float(frame_count - 1)
interp_pos = (1.0 - alpha) * reactant_pos + alpha * product_pos
⋮----
min_any = _min_distance_any_pair(frame_atoms)
⋮----
min_any_pair = min_any
min_any_pair_frame = frame_idx
⋮----
min_ads_surf = _min_adsorbate_surface_distance(frame_atoms, adsorbate_indices)
⋮----
min_ads_surface = min_ads_surf
min_ads_surface_frame = frame_idx
⋮----
min_any_pair = -1.0
⋮----
min_ads_surface = -1.0
⋮----
def _valid_indices(atoms: Atoms, indices: Iterable[int]) -> List[int]
⋮----
def _build_adsorbate_graph(atoms: Atoms, adsorbate_indices: Sequence[int]) -> Dict[int, set[int]]
⋮----
graph: Dict[int, set[int]] = {idx: set() for idx in adsorbate_indices}
⋮----
symbols = atoms.get_chemical_symbols()
⋮----
cutoff = _bond_cutoff(symbols[idx_i], symbols[idx_j])
dist = float(atoms.get_distance(idx_i, idx_j, mic=True))
⋮----
def _bond_cutoff(symbol_a: str, symbol_b: str) -> float
⋮----
radius_sum = covalent_radii[atomic_numbers[symbol_a]] + covalent_radii[atomic_numbers[symbol_b]]
⋮----
radius_sum = 1.6
⋮----
def _connected_components(indices: Sequence[int], graph: Dict[int, set[int]]) -> List[List[int]]
⋮----
remaining = set(indices)
components: List[List[int]] = []
⋮----
start = min(remaining)
queue = deque([start])
seen = {start}
⋮----
component = [start]
⋮----
current = queue.popleft()
⋮----
positions = np.array(atoms.get_positions(), dtype=float)
cell = np.array(atoms.cell, dtype=float)
⋮----
inv_t = np.linalg.pinv(cell.T)
⋮----
anchor_idx = _choose_component_anchor(component, symbols)
unwrapped: Dict[int, np.ndarray] = {anchor_idx: np.array(positions[anchor_idx], dtype=float)}
queue = deque([anchor_idx])
⋮----
current_frac = inv_t @ np.array(unwrapped[current], dtype=float)
⋮----
neighbor_frac = inv_t @ np.array(positions[neighbor], dtype=float)
delta = neighbor_frac - (inv_t @ np.array(positions[current], dtype=float))
⋮----
def _choose_component_anchor(component: Sequence[int], symbols: Sequence[str]) -> int
⋮----
non_h = [idx for idx in component if symbols[idx] != "H"]
⋮----
ref_frac = inv_t @ np.array(reference, dtype=float)
moving_frac = inv_t @ np.array(moving, dtype=float)
delta = ref_frac - moving_frac
shift = np.zeros(3, dtype=float)
⋮----
def _min_distance_any_pair(structure: Atoms) -> float
⋮----
dists = structure.get_all_distances(mic=False)
⋮----
value = float(np.min(dists))
⋮----
def _min_adsorbate_surface_distance(structure: Atoms, adsorbate_indices: Sequence[int]) -> float
⋮----
ads_set = set(adsorbate_indices)
surface_indices = [idx for idx in range(len(structure)) if idx not in ads_set]
⋮----
positions = structure.get_positions()
min_dist = float("inf")
⋮----
dist = float(np.linalg.norm(positions[ads_idx] - positions[surf_idx]))
min_dist = min(min_dist, dist)
````

## File: core/pathway/pathway_bundle.py
````python
"""Multi-pathway bundle builder and exporter.

Converts candidate pathways from the search engine into structured
output directories compatible with existing Agent4/5 consumption.
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
class PathwayBundleBuilder
⋮----
"""Build and export multi-pathway bundles to disk.

    Output layout::

        <output_base>/<run_id>/04_pathway/
            search_manifest.json
            iter_01_step_visualizations/
                path_00/
                    state_00_CO.vasp
                    step01_CO_to_CHO_reactant.vasp
                    step01_CO_to_CHO_product.vasp
                    ...
                path_01/
                    ...
    """
⋮----
"""Export a single pathway's structures and manifest entry.

        Args:
            pathway_id: unique pathway identifier
            pathway_index: numeric index (for path_XX directory)
            states: list of state dicts (species_label, structure_path, energy, etc.)
            steps: list of step dicts (reactant/product structure paths, step info)
            metadata: extra info to include in manifest

        Returns:
            Manifest entry dict for this pathway.
        """
path_dir = self._pathway_dir / f"path_{pathway_index:02d}"
⋮----
exported_states: List[Dict[str, Any]] = []
exported_steps: List[Dict[str, Any]] = []
⋮----
# Export state structures
⋮----
label = state.get("species_label", f"state_{idx}")
safe_label = _safe_name(label)
structure_path = state.get("structure_path")
⋮----
state_entry: Dict[str, Any] = {
⋮----
atoms = read(structure_path)
dest = path_dir / f"state_{idx:02d}_{safe_label}.vasp"
⋮----
# Export step structures (reactant/product pairs)
⋮----
step_name = step.get("step_name", f"step_{idx}")
safe_name = _safe_name(step_name)
⋮----
step_entry: Dict[str, Any] = {
⋮----
src_path = step.get(f"{endpoint}_structure_path")
⋮----
atoms = read(src_path)
dest = path_dir / f"step{idx + 1:02d}_{safe_name}_{endpoint}.vasp"
⋮----
manifest_entry = {
⋮----
# Write per-pathway manifest
manifest_file = path_dir / "pathway_manifest.json"
⋮----
"""Export multiple pathways and write a top-level search manifest.

        Args:
            pathways: list of dicts, each with:
                - pathway_id
                - is_retained
                - states (list)
                - steps (list)
                - metadata (optional)

        Returns:
            Complete search manifest dict.
        """
manifest_entries: List[Dict[str, Any]] = []
⋮----
entry = self.export_pathway(
⋮----
search_manifest = {
⋮----
# Write top-level manifest
manifest_path = self._pathway_dir.parent / "search_manifest.json"
⋮----
def _safe_name(name: str) -> str
⋮----
"""Convert a label to filesystem-safe name."""
⋮----
__all__ = ["PathwayBundleBuilder"]
````

## File: core/pathway/pathway_predictor.py
````python
"""
Pathway Predictor - 反应路径完整分析模块

这个模块整合了自由能预测和能垒预测，提供完整的反应路径分析功能。

对于给定的催化剂表面和反应路径，这个模块可以：
1. 预测每个反应中间体的吸附能（自由能）
2. 预测相邻中间体之间的反应能垒
3. 生成完整的能量剖面图
4. 识别速率决定步骤

使用方式:
    from pathway_predictor import PathwayPredictor

    # 初始化预测器
    predictor = PathwayPredictor(
        fairchem_root="/path/to/fairchem",
        model_name="uma-s-1p1",
        use_gpu=True,
    )

    # 方式1: 从吸附质列表预测完整反应路径
    result = predictor.predict_pathway(
        surface="Pt_111.vasp",
        adsorbates=["*O", "*OH", "*OOH", "*"],  # 按反应顺序
        calculate_barriers=True,  # 是否计算能垒
    )

    # 打印结果
    print(result.summary())

    # 绘制能量剖面图
    predictor.plot_energy_profile(result, output="energy_profile.png")

    # 方式2: 手动指定反应步骤
    result = predictor.predict_pathway_from_reactions(
        surface="Pt_111.vasp",
        reactions=[
            ("*O", "*OH", "O_hydrogenation"),
            ("*OH", "*OOH", "OH_hydrogenation"),
            ("*OOH", "*", "OOH_desorption"),
        ]
    )

作者: Claude
"""
⋮----
HAS_MATPLOTLIB = True
⋮----
HAS_MATPLOTLIB = False
⋮----
@dataclass
class PathwayStep
⋮----
"""反应路径中的单个步骤"""
step_index: int
name: str
reactant_adsorbate: str
product_adsorbate: str
reactant_energy: float  # eV
product_energy: float  # eV
reaction_energy: float  # eV
activation_energy: Optional[float] = None  # eV
transition_state_energy: Optional[float] = None  # eV
is_rate_determining: bool = False
⋮----
def __repr__(self)
⋮----
barrier_str = f", E_act={self.activation_energy:.3f} eV" if self.activation_energy else ""
⋮----
@dataclass
class CompletePathwayResult
⋮----
"""完整反应路径分析结果"""
surface_formula: str
adsorbates: List[str]
adsorbate_energies: Dict[str, float]  # adsorbate -> total energy (eV)
adsorption_energies: Dict[str, float]  # adsorbate -> E_ads (eV)
steps: List[PathwayStep]
rate_determining_step: Optional[PathwayStep]
max_barrier: Optional[float]  # eV
overall_reaction_energy: float  # eV
output_dir: str
⋮----
def summary(self) -> str
⋮----
"""生成结果摘要"""
lines = [
⋮----
e_ads = self.adsorption_energies.get(ads, 0.0)
e_total = self.adsorbate_energies.get(ads, 0.0)
⋮----
rds_marker = " [RDS]" if step.is_rate_determining else ""
⋮----
class PathwayPredictor
⋮----
"""
    反应路径完整分析预测器

    这个类整合了自由能预测和能垒预测，提供完整的反应路径分析。

    Parameters
    ----------
    fairchem_root : str
        Fairchem 项目根目录
    model_name : str, default="uma-s-1p1"
        使用的模型名称
    model_path : str, optional
        本地模型文件路径
    use_gpu : bool, default=True
        是否使用 GPU
    work_dir : str, optional
        工作目录
    keep_files : bool, default=True
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """
⋮----
# 原子参考能量 (eV) - 与 FairchemPredictor 保持一致
# 当原子从表面离开或加入时，使用这些能量进行校正
ATOMIC_REFERENCE_ENERGIES = {
⋮----
@staticmethod
    def parse_adsorbate_formula(formula: str) -> Counter
⋮----
"""
        解析吸附物化学式，返回元素组成

        支持格式: *CO, *CHOH, *CH3, *, CO, CHOH 等

        Parameters
        ----------
        formula : str
            吸附物化学式

        Returns
        -------
        Counter
            元素组成计数，如 {'C': 1, 'H': 1, 'O': 1} for *CHOH
        """
# 移除 * 前缀
formula = formula.strip().lstrip('*').strip()
⋮----
# 空字符串（纯表面）返回空计数
⋮----
# 解析化学式
# 匹配模式: 元素符号(大写+可选小写) + 可选数字
pattern = r'([A-Z][a-z]?)(\d*)'
composition = Counter()
⋮----
element = match.group(1)
count = int(match.group(2)) if match.group(2) else 1
if element:  # 忽略空匹配
⋮----
"""
        计算从反应物到产物的元素变化

        Parameters
        ----------
        reactant_formula : str
            反应物化学式 (如 *CHOH)
        product_formula : str
            产物化学式 (如 *CH)

        Returns
        -------
        Dict[str, int]
            元素变化量，正值表示产物多出的元素（需要从气相加入），
            负值表示反应物多出的元素（离开表面进入气相）
        """
reactant_comp = self.parse_adsorbate_formula(reactant_formula)
product_comp = self.parse_adsorbate_formula(product_formula)
⋮----
# 计算差异
change = {}
all_elements = set(reactant_comp.keys()) | set(product_comp.keys())
⋮----
diff = product_comp[elem] - reactant_comp[elem]
⋮----
"""
        计算气相能量校正

        当元素从表面离开时，需要加上该元素的参考能量（因为产物包含气相分子）
        当元素加入表面时，需要减去该元素的参考能量（因为反应物包含气相分子）

        对于反应 *CHOH → *CH:
        - 元素变化: O: -1, H: -1 (OH离开表面)
        - 校正 = -(-1) * E_O + -(-1) * E_H = E_O + E_H
        - 即产物侧需要加上离开的原子的能量

        Parameters
        ----------
        element_change : Dict[str, int]
            元素变化量 (product - reactant)

        Returns
        -------
        float
            气相能量校正 (eV)
        """
correction = 0.0
⋮----
# delta < 0 表示原子离开表面（产物有气相分子）
# delta > 0 表示原子加入表面（反应物有气相分子）
# 校正 = -delta * E_ref
# 这样：离开的原子会增加产物能量，加入的原子会增加反应物能量
⋮----
# 设置工作目录
⋮----
# 初始化日志
⋮----
# 延迟初始化子模块
⋮----
def _init_predictors(self)
⋮----
"""初始化子预测器"""
⋮----
model_path=self.model_path,  # 传递本地模型路径
⋮----
"""
        预测完整反应路径

        Parameters
        ----------
        surface : str or Atoms
            催化剂表面结构（可能包含吸附物）
        adsorbates : list of str or Atoms
            反应路径上的吸附质列表（按顺序）
        calculate_barriers : bool, default=True
            是否计算反应能垒
        num_sites : int, default=5
            每个吸附质尝试的位点数
        n_frames : int, default=10
            NEB 计算的帧数
        fmax : float, default=0.1
            力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        output_dir : str, optional
            输出目录
        current_adsorbate_indices : list of int, optional
            当前吸附物的原子索引（如果表面已有吸附物）

        Returns
        -------
        CompletePathwayResult
            完整反应路径分析结果
        """
⋮----
output_dir = os.path.join(self.work_dir, "pathway_complete")
⋮----
# 读取表面
⋮----
surface_atoms = read(surface)
⋮----
surface_atoms = surface.copy()
⋮----
surface_formula = surface_atoms.get_chemical_formula()
⋮----
# 如果表面已有吸附物，找到吸附位点并移除吸附物
adsorption_site = None
clean_surface = surface_atoms
⋮----
# 找到吸附物最下面的原子（z坐标最小）作为吸附位点
adsorbate_positions = surface_atoms.positions[current_adsorbate_indices]
lowest_atom_idx = np.argmin(adsorbate_positions[:, 2])
adsorption_site = adsorbate_positions[lowest_atom_idx]
⋮----
# 移除吸附物得到clean表面
clean_surface = surface_atoms.copy()
⋮----
# 重要：重新设置clean表面的tags，确保有正确的0/1分布
# 重构后的表面可能所有原子tags都是0（因为都被fixed了）
# 我们需要确保表面原子是tag=1，体相原子是tag=0
# 同时，FAIRChem要求原子顺序：bulk(0) -> surface(1) -> adsorbate(2)
⋮----
z_coords = clean_surface.positions[:, 2]
z_max = z_coords.max()
# 顶部3层原子作为表面（tag=1），其余为体相（tag=0）
is_surface = z_coords > z_max - 3.0
⋮----
# 关键修复：重新排序原子，使bulk原子在前，surface原子在后
# 这样tags才能正确对应到原子
bulk_indices = np.where(~is_surface)[0]
surface_indices = np.where(is_surface)[0]
⋮----
# 创建新的排序
new_order = np.concatenate([bulk_indices, surface_indices])
⋮----
# 重新排序：创建一个新的plain ASE Atoms对象
# 这避免了AtomsBatch等子类的内部状态问题
old_atoms = clean_surface
new_positions = old_atoms.positions[new_order]
new_symbols = [old_atoms.get_chemical_symbols()[i] for i in new_order]
⋮----
clean_surface = ASEAtoms(
⋮----
# 如果有constraints，也需要重新映射
⋮----
# 创建索引映射：old_index -> new_index
index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(new_order)}
new_constraints = []
⋮----
# 重新映射constraint的索引
⋮----
old_indices = constraint.index if hasattr(constraint.index, '__iter__') else [constraint.index]
new_indices = [index_map[i] for i in old_indices if i in index_map]
⋮----
# 设置tags：前N个是bulk(0)，后M个是surface(1)
n_bulk = len(bulk_indices)
n_surface = len(surface_indices)
new_tags = [0] * n_bulk + [1] * n_surface
⋮----
# 1. 预测所有吸附质的能量
⋮----
pathway_result = self._fairchem_predictor.predict_pathway_energies(
⋮----
fixed_adsorption_site=adsorption_site,  # 使用固定的吸附位点
⋮----
adsorbate_energies = pathway_result.adsorbate_energies
adsorption_energies = pathway_result.adsorption_energies
best_configs = pathway_result.best_configurations
⋮----
# 2. 计算相邻步骤之间的能垒
steps = []
previous_product_struct = None  # Track previous step's product for continuity
previous_product_adsorbate_indices = None
⋮----
reactant_ads = str(adsorbates[i])
product_ads = str(adsorbates[i + 1])
step_name = f"{reactant_ads}_to_{product_ads}"
⋮----
reactant_energy = adsorbate_energies[reactant_ads]
product_energy = adsorbate_energies[product_ads]
⋮----
# 计算元素变化并应用气相能量校正
element_change = self.calculate_element_change(reactant_ads, product_ads)
gas_correction = 0.0
⋮----
gas_correction = self.calculate_gas_phase_correction(element_change)
⋮----
# 反应能 = 产物能量 - 反应物能量 + 气相校正
# 气相校正已经包含了离开/加入表面的原子能量
reaction_energy = product_energy - reactant_energy + gas_correction
⋮----
# 获取反应物和产物结构
reactant_struct = best_configs[reactant_ads]
product_struct = best_configs[product_ads]
⋮----
# 计算吸附分子的原子索引（表面原子之后的所有原子）
# 使用clean_surface的长度，因为best_configs中的结构也是基于clean_surface的
surface_atom_count = len(clean_surface)
reactant_adsorbate_indices = list(range(surface_atom_count, len(reactant_struct)))
product_adsorbate_indices = list(range(surface_atom_count, len(product_struct)))
⋮----
# 计算能垒 - 传递前一步的产物信息以确保路径连续性
neb_result = self._barrier_predictor.predict_from_structures(
⋮----
relax_endpoints=False,  # 已经优化过了
⋮----
# Update previous product for next iteration
previous_product_struct = product_struct
previous_product_adsorbate_indices = product_adsorbate_indices
⋮----
activation_energy = neb_result.activation_energy_forward
# 只有当 transition_state_index 不为 None 时才获取 ts_energy
⋮----
ts_energy = neb_result.energies[neb_result.transition_state_index]
⋮----
ts_energy = None
⋮----
activation_energy = None
⋮----
# 创建步骤对象
step = PathwayStep(
⋮----
# 不计算能垒，只记录反应能
⋮----
# 3. 识别速率决定步骤（最大能垒）
rate_determining_step = None
max_barrier = None
⋮----
valid_barriers = [
⋮----
rate_determining_step = max(valid_barriers, key=lambda x: x[1])[0]
⋮----
max_barrier = rate_determining_step.activation_energy
⋮----
# 4. 计算总反应能
overall_reaction_energy = (
⋮----
"""
        从预计算的反应物/产物结构计算能垒

        适用于已有 Agent4/5 生成的端点结构，跳过吸附能计算阶段，
        直接对每步运行 NEB 并生成完整路径结果。

        Parameters
        ----------
        steps : list of dict
            每步包含:
              - name: str, 步骤名称
              - reactant: Atoms, 反应物结构
              - product: Atoms, 产物结构
              - reactant_formula: str, 反应物标签 (如 "*CO")
              - product_formula: str, 产物标签 (如 "*CHO")
        n_frames : int, default=10
            NEB 帧数
        fmax : float, default=0.1
            CI-NEB 力收敛标准 (eV/Å)
        max_steps : int, default=300
            最大优化步数
        relax_endpoints : bool, default=False
            是否先优化端点结构
        output_dir : str, optional
            输出目录

        Returns
        -------
        CompletePathwayResult
            完整反应路径分析结果
        """
⋮----
output_dir = os.path.join(self.work_dir, "pathway_from_structures")
⋮----
# Collect adsorbate labels
adsorbates = [steps[0]["reactant_formula"]]
⋮----
# Surface formula from first reactant
surface_formula = steps[0]["reactant"].get_chemical_formula()
⋮----
# Placeholder energies (NEB gives relative, not absolute)
adsorbate_energies = {}
adsorption_energies = {}
⋮----
pathway_steps = []
previous_product_struct = None
previous_product_ads_indices = None
⋮----
name = step["name"]
reactant = step["reactant"]
product = step["product"]
reactant_label = step["reactant_formula"]
product_label = step["product_formula"]
⋮----
# Detect adsorbate indices (non-surface atoms)
reactant_ads = list(step.get("reactant_adsorbate_indices", []))
product_ads = list(step.get("product_adsorbate_indices", []))
⋮----
# Staged (parked) atom indices for constrained endpoint relaxation
r_staged = list(step.get("reactant_staged_indices", []))
p_staged = list(step.get("product_staged_indices", []))
⋮----
# Gas-phase correction for element changes
element_change = self.calculate_element_change(reactant_label, product_label)
⋮----
reaction_energy = 0.0
⋮----
previous_product_struct = product
previous_product_ads_indices = product_ads
⋮----
reaction_energy = neb_result.reaction_energy
⋮----
# Store energies for profile
reactant_e = neb_result.energies[0] if activation_energy is not None else 0.0
product_e = neb_result.energies[-1] if activation_energy is not None else 0.0
⋮----
# Identify rate-determining step
⋮----
# Overall reaction energy
first_e = adsorbate_energies.get(adsorbates[0], 0.0)
last_e = adsorbate_energies.get(adsorbates[-1], 0.0)
overall_reaction_energy = last_e - first_e
⋮----
"""
        从反应列表预测路径

        Parameters
        ----------
        surface : str or Atoms
            催化剂表面
        reactions : list of tuple
            反应列表，每个元组为 (reactant, product, name)
        num_sites : int, default=5
            位点数
        n_frames : int, default=10
            NEB 帧数
        fmax : float, default=0.1
            力收敛标准
        max_steps : int, default=300
            最大优化步数
        output_dir : str, optional
            输出目录

        Returns
        -------
        CompletePathwayResult
            完整反应路径分析结果
        """
# 提取所有唯一的吸附质
unique_adsorbates = set()
⋮----
# 先按反应顺序排列吸附质
adsorbates_ordered = []
⋮----
# 使用标准路径预测
⋮----
"""
        绘制能量剖面图（使用 EnergyDiagramPlotter）

        从每步的 activation_energy 和 reaction_energy 构建累积能量剖面，
        而非使用绝对能量（不同步骤原子数不同，绝对能量不可比）。

        Parameters
        ----------
        result : CompletePathwayResult
            路径分析结果
        output : str, optional
            输出文件路径
        figsize : tuple, default=(12, 6)
            图像大小
        show_barriers : bool, default=True
            是否显示能垒
        """
⋮----
# Build cumulative energy profile from relative barriers
energies = [0.0]
labels = [result.adsorbates[0]]
is_ts = [False]
⋮----
cumulative_energy = 0.0
⋮----
ea_fwd = step.activation_energy if step.activation_energy is not None else 0.0
reaction_e = step.reaction_energy if step.reaction_energy is not None else 0.0
⋮----
# TS energy = current level + forward barrier
ts_energy = cumulative_energy + ea_fwd
⋮----
# Product energy = current level + reaction energy
⋮----
style = DiagramStyle()
⋮----
plotter = EnergyDiagramPlotter(style=style)
⋮----
def __del__(self)
⋮----
"""清理临时文件"""
````

## File: core/pathway/pruning_policy.py
````python
"""Pruning policy for mechanism search beam management.

Implements sibling comparison, beam selection, and budget control.
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
@dataclass
class PruningPolicy
⋮----
"""Configurable policy for pathway branch pruning.

    The core logic:
    - If sibling energy gap <= ``delta_keep`` → keep both
    - If sibling energy gap >= ``delta_prune`` → prune the higher one
    - In between → keep but mark as marginal
    - After sibling comparison, apply beam width to limit frontier size
    """
⋮----
delta_keep: float = 0.10   # eV
delta_prune: float = 0.30  # eV
beam_width: int = 8
max_depth: int = 8
max_state_evaluations: int = 40
⋮----
"""Compare sibling candidates and decide keep/prune.

        Args:
            sibling_estimates: list of dicts, each with at least
                ``state_id`` and ``free_energy_eV``.

        Returns:
            Dict with ``kept``, ``pruned``, ``marginal`` lists of state_ids,
            plus ``decisions`` mapping state_id → reason.
        """
⋮----
# Sort by energy (None → infinity)
def _energy(e: Dict) -> float
⋮----
v = e.get("free_energy_eV")
⋮----
sorted_siblings = sorted(sibling_estimates, key=_energy)
best_energy = _energy(sorted_siblings[0])
⋮----
kept: List[str] = []
pruned: List[str] = []
marginal: List[str] = []
decisions: Dict[str, str] = {}
⋮----
sid = est.get("state_id", "?")
energy = _energy(est)
gap = energy - best_energy
⋮----
"""Select the top ``beam_width`` nodes from a frontier.

        Args:
            nodes: list of dicts with at least ``state_id`` and ``free_energy_eV``.

        Returns:
            (selected, dropped) tuple.
        """
def _energy(n: Dict) -> float
⋮----
v = n.get("free_energy_eV")
⋮----
sorted_nodes = sorted(nodes, key=_energy)
selected = sorted_nodes[: self.beam_width]
dropped = sorted_nodes[self.beam_width :]
⋮----
"""Check whether the search should terminate."""
⋮----
__all__ = ["PruningPolicy"]
````

## File: core/pathway/search_engine.py
````python
"""Mechanism search engine: best-first + beam search over catalytic states.

Expands candidate states layer by layer, evaluates free energies, prunes
siblings, and extracts complete pathways from root to goal.
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
# ---------------------------------------------------------------------------
# Search node
⋮----
@dataclass
class PathSearchNode
⋮----
"""One node in the search tree."""
state_id: str
species_label: str
elements: Dict[str, int] = field(default_factory=dict)
depth: int = 0
parent_id: Optional[str] = None
step_from_parent: Optional[Dict[str, Any]] = None
free_energy_eV: Optional[float] = None
is_goal: bool = False
is_pruned: bool = False
prune_reason: str = ""
⋮----
# Search engine
⋮----
class MechanismSearchEngine
⋮----
"""Layer-by-layer beam search over catalytic state space."""
⋮----
# Search state
⋮----
self._visited_labels: Set[str] = set()  # species_label dedup across depths
⋮----
# ----- public API -----
⋮----
"""Create the root node for the search."""
⋮----
elements = parse_species_elements(species_label)
⋮----
state_id = self._make_state_id(species_label, depth=0)
root = PathSearchNode(
⋮----
def expand_node(self, node_id: str) -> List[str]
⋮----
"""Expand a single node into its children."""
node = self._nodes.get(node_id)
⋮----
candidates = self.generator.generate_candidates(
⋮----
child_ids: List[str] = []
⋮----
product_label = cand.get("product_label", "")
product_elements = cand.get("product_elements", {})
⋮----
product_elements = parse_species_elements(product_label)
⋮----
child_id = self._make_state_id(product_label, depth=node.depth + 1)
# Skip duplicates (by state_id AND by normalized label across depths)
normalized = self._strip_phase(product_label)
⋮----
child = PathSearchNode(
⋮----
"""Evaluate free energies for all frontier nodes."""
count = 0
⋮----
node = self._nodes.get(nid)
⋮----
# Skip gas-phase and composite labels (UMA can't handle them)
label = node.species_label
⋮----
node.free_energy_eV = 0.0  # placeholder for gas/composite
⋮----
extra_kwargs: Dict[str, Any] = {}
⋮----
est = self.energy_router.estimate_state_free_energy(
⋮----
# Release GPU memory after each frontier evaluation batch
⋮----
def prune_frontier(self) -> Dict[str, Any]
⋮----
"""Apply sibling comparison + beam width to the current frontier."""
# Group siblings by parent
parent_groups: Dict[Optional[str], List[str]] = {}
⋮----
node = self._nodes[nid]
parent = node.parent_id
⋮----
all_kept: List[str] = []
all_pruned: List[str] = []
decisions: Dict[str, str] = {}
⋮----
estimates = []
⋮----
node = self._nodes[cid]
⋮----
result = self.policy.compare_siblings(estimates)
kept = result["kept"] + result["marginal"]
pruned = result["pruned"]
⋮----
pnode = self._nodes.get(pid)
⋮----
# Apply beam width
kept_nodes = [{"state_id": k, "free_energy_eV": self._nodes[k].free_energy_eV} for k in all_kept]
⋮----
selected_ids = [s["state_id"] for s in selected]
⋮----
did = d["state_id"]
dnode = self._nodes.get(did)
⋮----
"""Check if a node matches the target state.

        Handles phase-agnostic matching: ``*CH4`` matches ``CH4(g)`` because
        they share the same element composition.
        """
⋮----
# Label match (normalized — strips *, (g), (s), whitespace)
⋮----
# Element match (phase-agnostic)
⋮----
# Phase-stripped label match: *CH4 vs CH4(g) → both normalize to "ch4"
⋮----
def backtrack_path(self, goal_id: str) -> List[PathSearchNode]
⋮----
"""Trace back from a goal node to the root, returning the path."""
path: List[PathSearchNode] = []
current_id: Optional[str] = goal_id
⋮----
node = self._nodes.get(current_id)
⋮----
current_id = node.parent_id
⋮----
"""Run the full search loop.

        Returns a dict with ``pathways`` (list of node sequences),
        ``stats``, and ``all_nodes``.
        """
⋮----
target_elements = parse_species_elements(target_label)
⋮----
depth = 0
⋮----
# Check stopping condition
⋮----
# Expand all frontier nodes
new_children: List[str] = []
⋮----
children = self.expand_node(nid)
⋮----
# Check for goal states among children
⋮----
# Move frontier to children
⋮----
# Evaluate and prune
⋮----
# Extract pathways from goal nodes
pathways = []
⋮----
path = self.backtrack_path(gid)
⋮----
# If no goal reached, extract best frontier paths as partial candidates
⋮----
frontier_nodes = [
# Sort by energy (lowest first)
⋮----
# Take top-3 best frontier paths
⋮----
path = self.backtrack_path(nid)
⋮----
# ----- graph enumeration for systematic search -----
⋮----
"""BFS enumeration of ALL paths from root to target (no energy eval).

        Unlike ``search_until_stop``, this does NOT prune by energy.  It builds
        a full reachability graph and extracts all distinct paths to the target,
        suitable for subsequent LLM filtering + UMA evaluation.

        Each returned path is a list of dicts: ``{"species_label", "elements", "step_info"}``.
        """
⋮----
# Build adjacency graph: label → [(child_label, child_elements, step_info)]
graph: Dict[str, List[Dict[str, Any]]] = {}
expand_queue: List[Tuple[str, Dict[str, int], int]] = []
⋮----
# Seed from root
⋮----
root_node = self._nodes[self._frontier[0]]
⋮----
visited_for_expand: Set[str] = {self._strip_phase(root_node.species_label)}
⋮----
children = self.generator.generate_candidates(
adj = []
⋮----
prod = c.get("product_label", "")
prod_elem = c.get("product_elements", {})
⋮----
prod_elem = parse_species_elements(prod)
⋮----
norm = self._strip_phase(prod)
⋮----
# DFS to find all paths from root to target
root_label = root_node.species_label
all_paths: List[List[Dict[str, Any]]] = []
⋮----
def _dfs(current: str, current_elem: Dict[str, int], path: List[Dict[str, Any]], visited: Set[str])
⋮----
# Check goal
⋮----
child_label = child["label"]
child_norm = self._strip_phase(child_label)
⋮----
continue  # no cycles within same path
⋮----
root_entry = {
⋮----
# ----- internal helpers -----
⋮----
@staticmethod
    def _make_state_id(label: str, depth: int) -> str
⋮----
"""Create a unique state ID."""
safe = label.replace("*", "ads_").replace("(", "_").replace(")", "_").replace("+", "_")
⋮----
@staticmethod
    def _normalize_label(label: str) -> str
⋮----
@staticmethod
    def _strip_phase(label: str) -> str
⋮----
"""Remove phase markers: '*CO' → 'co', 'CH4(g)' → 'ch4'."""
s = label.strip().lower().replace(" ", "")
s = s.replace("*", "").replace("(g)", "").replace("(s)", "").replace("(l)", "")
⋮----
__all__ = ["PathSearchNode", "MechanismSearchEngine"]
````

## File: core/reconstruction/__init__.py
````python

````

## File: core/reconstruction/adsorbdiff_predictor.py
````python
"""
AdsorbDiff Predictor - 吸附位点预测封装类

这个模块提供了一个完整的封装类，用于预测最优吸附位点和取向。

使用方式:
    from adsorbdiff_predictor import AdsorbDiffPredictor

    # 初始化预测器
    predictor = AdsorbDiffPredictor(
        adsorbdiff_root="/path/to/AdsorbDiff",
        checkpoint_path="/path/to/checkpoint.pt",  # 可选，默认使用PaiNN模型
        use_gpu=True,
    )

    # 方式1: 表面+吸附物分别提供
    result = predictor.predict(
        surface="slab.vasp",
        adsorbate="*CO",  # 可以是SMILES字符串或ASE Atoms
        has_adsorbate=False,
    )

    # 方式2: 已有吸附物的表面结构
    result = predictor.predict(
        surface="ads_slab.vasp",
        has_adsorbate=True,  # 吸附物已在结构中
        adsorbate_tag=2,  # 吸附物原子的tag
    )

作者: Claude
"""
⋮----
@dataclass
class AdsorptionSite
⋮----
"""单个吸附位点信息"""
site_id: int
position: np.ndarray  # 吸附位点坐标
binding_type: str  # "ontop", "bridge", "hollow", "random"
energy: Optional[float] = None  # 预测能量 (eV)
⋮----
@dataclass
class AdsorptionResult
⋮----
"""吸附预测结果"""
initial_structure: "ase.Atoms"  # 初始结构
final_structure: "ase.Atoms"  # 最终弛豫结构
adsorbate_position: np.ndarray  # 最终吸附位置 (COM)
adsorbate_orientation: np.ndarray  # 最终取向
energy: float  # 预测能量 (eV)
forces_max: float  # 最大残余力 (eV/Å)
trajectory_path: Optional[str] = None  # 轨迹文件路径
site_info: Optional[AdsorptionSite] = None  # 初始位点信息
⋮----
def __repr__(self)
⋮----
@dataclass
class PredictionOutput
⋮----
"""完整的预测输出"""
surface_formula: str
adsorbate_formula: str
num_sites_sampled: int
results: List[AdsorptionResult]
best_result: AdsorptionResult
output_dir: Optional[str] = None
use_proxy_energy: bool = False  # 是否使用了代理能量
⋮----
def get_top_n(self, n: int) -> List[AdsorptionResult]
⋮----
"""获取能量最低的前N个结果"""
sorted_results = sorted(self.results, key=lambda x: x.energy)
⋮----
def summary(self, top_n: int = 5) -> str
⋮----
"""生成预测结果摘要"""
energy_label = "Proxy Energy" if self.use_proxy_energy else "Energy"
lines = [
⋮----
class AdsorbDiffPredictor
⋮----
"""
    AdsorbDiff吸附位点预测器

    这个类封装了AdsorbDiff项目的所有功能，提供简洁的API用于预测最优吸附位点和取向。

    Parameters
    ----------
    adsorbdiff_root : str
        AdsorbDiff项目的根目录路径
    checkpoint_path : str, optional
        扩散模型checkpoint文件路径。如果不指定，将使用默认的PaiNN模型
    relax_checkpoint_path : str, optional
        弛豫模型checkpoint文件路径（用于能量和力的计算）。
        如果不指定，将使用代理能量进行配置排序（基于吸附位置而非实际能量计算）
    use_gpu : bool, default=True
        是否使用GPU进行推理
    cutoff : float, default=12.0
        截断半径 (Å)
    max_neighbors : int, default=50
        最大邻居数
    num_sites : int, default=10
        采样的吸附位点数量
    num_augmentations_per_site : int, default=1
        每个位点的增强数量
    placement_mode : str, default="heuristic"
        吸附物放置模式: "random", "heuristic", "random_site_heuristic_placement"
    interstitial_gap : float, default=2.0
        吸附物与表面的最小距离 (Å)
    diffusion_steps : int, default=100
        扩散步数
    ads_std_low : float, default=0.1
        吸附物位置噪声的最小标准差
    ads_std_high : float, default=10.0
        吸附物位置噪声的最大标准差
    rot_std_low : float, default=0.01
        旋转噪声的最小标准差
    rot_std_high : float, default=1.55
        旋转噪声的最大标准差
    relax_steps : int, default=100
        弛豫优化最大步数
    relax_fmax : float, default=0.05
        弛豫收敛力阈值 (eV/Å)
    seed : int, optional
        随机种子
    work_dir : str, optional
        工作目录。如果不指定，将创建临时目录
    keep_files : bool, default=False
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """
⋮----
interstitial_gap: float = 0.5,  # Gap between adsorbate and surface (was 2.0, too high)
⋮----
# 设置checkpoint路径
⋮----
default_painn = os.path.join(self.adsorbdiff_root, "PT_zeroshot_painn.pt")
default_eqv2 = os.path.join(self.adsorbdiff_root, "PT_fewshot_eqv2_cond.pt")
⋮----
# 设置工作目录
⋮----
# 添加AdsorbDiff路径到系统路径
⋮----
# 设置随机种子
⋮----
# 验证环境
⋮----
# 延迟加载模型（第一次预测时加载）
⋮----
def _validate_environment(self)
⋮----
"""验证运行环境"""
⋮----
# 检查必要的数据库文件
ads_db = os.path.join(self.adsorbdiff_root, "adsorbdiff/placement/pkls/adsorbates.pkl")
⋮----
def _log(self, message: str)
⋮----
"""输出日志"""
⋮----
def _load_diffusion_calculator(self)
⋮----
"""加载扩散模型"""
⋮----
def _load_relax_calculator(self)
⋮----
"""加载弛豫模型"""
⋮----
# 尝试使用fairchem的OCPCalculator
⋮----
# 回退到AdsorbDiff的calculator
⋮----
# 没有弛豫模型时返回None，使用代理能量
⋮----
def _convert_surface(self, surface) -> "ase.Atoms"
⋮----
"""将输入表面转换为ASE Atoms"""
⋮----
# 尝试从pymatgen Structure转换
⋮----
def _convert_adsorbate(self, adsorbate) -> "ase.Atoms"
⋮----
"""将输入吸附物转换为ASE Atoms"""
⋮----
# 尝试作为SMILES解析
⋮----
# 从数据库查找
⋮----
# 从数据库按ID查找
⋮----
def _get_adsorbate_from_smiles(self, smiles: str) -> "ase.Atoms"
⋮----
"""从数据库获取吸附物（按SMILES）"""
ads_db_path = os.path.join(
⋮----
ads_db = pickle.load(f)
⋮----
# 搜索匹配的SMILES
⋮----
# 如果没找到，列出可用的吸附物
available = [ads_info[1] for ads_info in ads_db.values()][:20]
⋮----
def _get_adsorbate_from_db(self, idx: int) -> "ase.Atoms"
⋮----
"""从数据库获取吸附物（按ID）"""
⋮----
"""从表面结构中提取吸附物"""
tags = atoms.get_tags()
ads_mask = tags == adsorbate_tag
slab_mask = ~ads_mask
⋮----
slab = atoms[slab_mask]
adsorbate = atoms[ads_mask]
⋮----
# 设置正确的tags
⋮----
def _tag_surface_atoms(self, atoms: "ase.Atoms") -> "ase.Atoms"
⋮----
"""为表面原子添加tags (1=表面, 0=subsurface)"""
⋮----
# 已经有tags
⋮----
# 使用简单的高度判断
positions = atoms.get_positions()
z_coords = positions[:, 2]
z_max = z_coords.max()
z_threshold = z_max - 2.5  # 表面2.5 Å以内的原子
⋮----
new_tags = np.where(z_coords > z_threshold, 1, 0)
⋮----
"""创建吸附物-表面配置"""
⋮----
# 创建临时Bulk对象
class DummyBulk
⋮----
def __init__(self, atoms)
⋮----
bulk = DummyBulk(slab_atoms)
⋮----
# 创建Slab对象（绕过验证）
slab = object.__new__(Slab)
⋮----
# 确保slab有surface tags
⋮----
# 创建Adsorbate对象
adsorbate = Adsorbate(
⋮----
# 创建吸附配置
config = AdsorbateSlabConfig(
⋮----
"""运行扩散采样"""
calc = self._load_diffusion_calculator()
⋮----
"""运行结构弛豫"""
# 如果没有弛豫模型，使用代理能量
⋮----
calc = self._load_relax_calculator()
⋮----
opt = BFGS(adslab, trajectory=traj_path, logfile=None)
⋮----
opt = BFGS(adslab, logfile=None)
⋮----
energy = adslab.get_potential_energy()
forces = adslab.get_forces()
fmax = np.max(np.linalg.norm(forces, axis=1))
⋮----
"""
        计算代理能量（当没有弛豫模型时使用）

        使用吸附物与表面的距离和位置作为代理能量指标。
        较低的能量表示更稳定的配置。
        """
tags = adslab.get_tags()
ads_mask = tags == 2
surf_mask = tags == 1
⋮----
ads_positions = adslab.positions[ads_mask]
surf_positions = adslab.positions[surf_mask]
⋮----
# 计算吸附物质心
ads_com = ads_positions.mean(axis=0)
⋮----
# 计算到表面原子的最小距离
distances = []
⋮----
dist = np.linalg.norm(ads_com[:2] - surf_pos[:2])  # xy平面距离
⋮----
min_xy_dist = min(distances)
⋮----
# 计算吸附物到表面的高度
surf_z_max = surf_positions[:, 2].max()
ads_z_min = ads_positions[:, 2].min()
height = ads_z_min - surf_z_max
⋮----
# 代理能量：距离表面原子越近、高度在合理范围内越好
# 理想高度约1.5-2.5 Å
ideal_height = 2.0
height_penalty = (height - ideal_height) ** 2
⋮----
# 代理能量 (越低越好)
proxy_energy = min_xy_dist * 0.5 + height_penalty * 0.3
⋮----
# 代理力 (使用一个固定的小值)
proxy_fmax = 0.1
⋮----
def _get_adsorbate_info(self, atoms: "ase.Atoms") -> Tuple[np.ndarray, np.ndarray]
⋮----
"""获取吸附物的位置和取向信息"""
⋮----
ads_positions = atoms.positions[ads_mask]
⋮----
# 计算质心
com = ads_positions.mean(axis=0)
⋮----
# 计算主轴方向（简化）
⋮----
# 使用第一个到最后一个原子的方向
orientation = ads_positions[-1] - ads_positions[0]
norm = np.linalg.norm(orientation)
⋮----
orientation = orientation / norm
⋮----
orientation = np.array([0, 0, 1])
⋮----
"""
        预测最优吸附位点和取向

        Parameters
        ----------
        surface : str or ase.Atoms or pymatgen.Structure
            表面结构。可以是文件路径、ASE Atoms对象或pymatgen Structure
        adsorbate : str or int or ase.Atoms, optional
            吸附物。可以是:
            - SMILES字符串 (如 "*CO", "*O", "*OH")
            - 数据库ID (整数)
            - ASE Atoms对象
            - 文件路径
            如果has_adsorbate=True，则忽略此参数
        has_adsorbate : bool, default=False
            表面结构中是否已包含吸附物
        adsorbate_tag : int, default=2
            当has_adsorbate=True时，用于识别吸附物原子的tag值
        binding_indices : list of int, optional
            吸附物的结合原子索引列表
        num_samples : int, optional
            采样的配置数量。如果不指定，使用初始化时的num_sites
        output_dir : str, optional
            输出目录
        save_trajectory : bool, default=True
            是否保存扩散轨迹

        Returns
        -------
        PredictionOutput
            预测结果，包含最优配置和所有采样配置
        """
⋮----
# 设置输出目录
⋮----
output_dir = os.path.join(self.work_dir, "prediction")
⋮----
# 转换表面结构
surface_atoms = self._convert_surface(surface)
⋮----
# 处理吸附物
⋮----
# 从表面提取吸附物
⋮----
slab_atoms = surface_atoms
adsorbate_atoms = self._convert_adsorbate(adsorbate)
⋮----
# 确保表面有正确的tags
slab_atoms = self._tag_surface_atoms(slab_atoms)
⋮----
# 使用已有的吸附物位置作为初始配置
adslab_list = [surface_atoms.copy()]
metadata_list = [{"site": self._get_adsorbate_info(surface_atoms)[0]}]
⋮----
# 运行扩散和弛豫
results = []
traj_dir = os.path.join(output_dir, "trajectories")
⋮----
initial_atoms = adslab.copy()
⋮----
# 运行扩散
diff_traj_path = os.path.join(traj_dir, f"diff_{i}")
⋮----
diffused_adslab = self._run_diffusion(adslab, diff_traj_path)
⋮----
diffused_adslab = adslab.copy()
⋮----
# 获取扩散后的吸附位置
⋮----
# 运行弛豫
relax_traj_path = os.path.join(traj_dir, f"relax_{i}.traj") if save_trajectory else None
⋮----
# 获取最终位置
⋮----
result = AdsorptionResult(
⋮----
# 找到最优配置
best_result = min(results, key=lambda x: x.energy)
⋮----
output = PredictionOutput(
⋮----
# 保存最优结构
⋮----
best_path = os.path.join(output_dir, "best_config.vasp")
⋮----
"""
        对单个吸附物-表面配置进行预测

        Parameters
        ----------
        adslab : ase.Atoms
            已经放置好吸附物的表面结构（吸附物tag=2）
        output_dir : str, optional
            输出目录

        Returns
        -------
        AdsorptionResult
            预测结果
        """
⋮----
output_dir = os.path.join(self.work_dir, "single_prediction")
⋮----
traj_dir = os.path.join(output_dir, "trajectory")
⋮----
diffused_adslab = self._run_diffusion(adslab, traj_dir)
⋮----
relax_traj_path = os.path.join(output_dir, "relax.traj")
⋮----
def list_available_adsorbates(self) -> List[Dict]
⋮----
"""列出数据库中所有可用的吸附物"""
⋮----
adsorbates = []
⋮----
def __del__(self)
⋮----
"""清理临时目录"""
⋮----
# 便捷函数
⋮----
"""
    便捷函数：预测最优吸附位点

    Parameters
    ----------
    surface : str or ase.Atoms
        表面结构
    adsorbate : str or int or ase.Atoms
        吸附物 (SMILES, 数据库ID, 或ASE Atoms)
    adsorbdiff_root : str
        AdsorbDiff项目根目录
    checkpoint_path : str, optional
        模型checkpoint路径
    use_gpu : bool, default=True
        是否使用GPU
    num_samples : int, default=10
        采样配置数量
    output_dir : str, optional
        输出目录
    **kwargs
        其他参数传递给AdsorbDiffPredictor

    Returns
    -------
    PredictionOutput
        预测结果
    """
predictor = AdsorbDiffPredictor(
⋮----
parser = argparse.ArgumentParser(description="AdsorbDiff Adsorption Site Predictor")
⋮----
args = parser.parse_args()
⋮----
result = predictor.predict(
````

## File: core/reconstruction/cp_mace_predictor.py
````python
"""
CP-MACE Predictor - 常电位分子动力学模拟封装类

这个模块提供了一个完整的封装类，用于训练CP-MACE模型并进行常电位分子动力学模拟。

CP-MACE (Constant-Potential MACE) 是MACE框架的扩展，支持在巨正则系综条件下
对电化学界面进行常电位分子模拟。

使用方式:
    from cp_mace_predictor import CPMACEPredictor

    # 初始化预测器
    predictor = CPMACEPredictor(
        cp_mace_root="/path/to/CP-MACE",
        use_gpu=True,
    )

    # 方式1: 使用预训练模型直接模拟
    result = predictor.simulate(
        structure="init.xyz",  # 带有吸附物的表面结构
        model_paths=["model1.model", "model2.model"],
        target_potential=-3.36,  # 目标电极电位 (V)
        temperature=300.0,
        steps=1000,
    )

    # 方式2: 先训练模型，再模拟
    train_result = predictor.train(
        train_file="train.xyz",  # 训练数据文件
        model_name="my_model",
        max_num_epochs=300,
    )

    sim_result = predictor.simulate(
        structure="init.xyz",
        model_paths=[train_result.model_path],
        target_potential=-3.36,
    )

作者: Claude
"""
⋮----
@dataclass
class TrainingResult
⋮----
"""训练结果"""
model_name: str
model_path: str
checkpoint_dir: str
train_file: str
valid_file: Optional[str]
num_epochs: int
final_loss: Optional[float]
log_file: str
⋮----
def __repr__(self)
⋮----
@dataclass
class SimulationStep
⋮----
"""单步模拟信息"""
step: int
time_fs: float
energy: float  # eV
kinetic_energy: float  # eV
potential_energy: float  # eV
temperature: float  # K
fermi_level: float  # V
force_std: float  # eV/Å
fermi_std: float  # V
electron_number: float
⋮----
@dataclass
class SimulationResult
⋮----
"""模拟结果"""
surface_formula: str
model_paths: List[str]
target_potential: float
temperature: float
total_steps: int
completed_steps: int
timestep_fs: float
trajectory: List["ase.Atoms"]
energy_history: List[float]
fermi_history: List[float]
electron_history: List[float]
force_std_history: List[float]
fermi_std_history: List[float]
final_structure: "ase.Atoms"
output_dir: str
high_uncertainty_structures: List["ase.Atoms"] = field(default_factory=list)
⋮----
def summary(self, show_last_n: int = 10) -> str
⋮----
"""生成模拟结果摘要"""
lines = [
⋮----
class CPMACEPredictor
⋮----
"""
    CP-MACE 常电位分子动力学预测器

    这个类封装了CP-MACE项目的所有功能，提供简洁的API用于训练CP-MACE模型
    并进行常电位分子动力学模拟。

    Parameters
    ----------
    cp_mace_root : str
        CP-MACE项目的根目录路径
    use_gpu : bool, default=True
        是否使用GPU进行训练和推理
    device : str, optional
        指定设备 ("cuda" 或 "cpu")，默认根据use_gpu自动选择
    default_dtype : str, default="float64"
        默认数据类型
    work_dir : str, optional
        工作目录。如果不指定，将创建临时目录
    keep_files : bool, default=True
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """
⋮----
# 设置工作目录
⋮----
# 添加CP-MACE路径到系统路径
mace_path = os.path.join(self.cp_mace_root, "mace")
⋮----
simulation_path = os.path.join(self.cp_mace_root, "simulation")
⋮----
# 验证环境
⋮----
# 延迟加载
⋮----
def _validate_environment(self)
⋮----
"""验证运行环境"""
⋮----
mace_init = os.path.join(self.cp_mace_root, "mace", "__init__.py")
⋮----
simulation_dir = os.path.join(self.cp_mace_root, "simulation")
⋮----
def _log(self, message: str)
⋮----
"""输出日志"""
⋮----
def _convert_structure(self, structure) -> "ase.Atoms"
⋮----
"""将输入结构转换为ASE Atoms"""
⋮----
# 尝试从pymatgen Structure转换
⋮----
"""
        准备CP-MACE训练数据

        CP-MACE需要xyz格式的训练数据，每个结构需要包含electron和potential标签。

        Parameters
        ----------
        structures : list
            结构列表 (ASE Atoms或文件路径)
        electron_numbers : list of float
            每个结构的电子数或净电荷
        potentials : list of float
            每个结构的费米能级/电极电位
        energies : list of float
            每个结构的参考能量 (eV)
        forces : list of np.ndarray
            每个结构的参考原子力 (eV/Å)
        output_file : str
            输出xyz文件路径

        Returns
        -------
        str
            输出文件路径
        """
⋮----
atoms_list = []
⋮----
atoms = self._convert_structure(struct)
⋮----
# 设置参考数据
⋮----
# 写入extended xyz格式
output_path = os.path.abspath(output_file)
⋮----
"""
        训练CP-MACE模型

        Parameters
        ----------
        train_file : str
            训练数据文件路径 (extended xyz格式)
        model_name : str, default="CP_MACE_model"
            模型名称
        valid_file : str, optional
            验证数据文件路径。如果不指定，使用valid_fraction从训练数据中划分
        valid_fraction : float, default=0.05
            验证集比例 (当valid_file未指定时使用)
        model_type : str, default="FermiMACE"
            模型类型: "FermiMACE" (推荐，节点增强方法) 或 "FermiMACE_2" (全局状态方法)
        hidden_irreps : str, default="128x0e + 128x1o"
            隐藏层不可约表示
        r_max : float, default=5.0
            截断半径 (Å)
        batch_size : int, default=10
            批次大小
        max_num_epochs : int, default=300
            最大训练轮数
        energy_weight : float, default=1.0
            能量损失权重
        forces_weight : float, default=100.0
            力损失权重
        potential_weight : float, default=10.0
            电位/费米能级损失权重
        E0s : str, default="average"
            E0计算方式
        ema : bool, default=True
            是否使用指数移动平均
        ema_decay : float, default=0.99
            EMA衰减率
        amsgrad : bool, default=True
            是否使用AMSGrad优化器
        seed : int, default=1
            随机种子
        restart_latest : bool, default=True
            是否从最新检查点重启
        output_dir : str, optional
            输出目录
        extra_args : dict, optional
            额外的命令行参数

        Returns
        -------
        TrainingResult
            训练结果
        """
⋮----
# 设置输出目录
⋮----
output_dir = os.path.join(self.work_dir, "training", model_name)
⋮----
train_file = os.path.abspath(train_file)
⋮----
# 构建命令
cmd = [
⋮----
# 添加额外参数
⋮----
# 运行训练
log_file = os.path.join(output_dir, "train.log")
⋮----
result = subprocess.run(
⋮----
# 查找生成的模型文件
model_path = os.path.join(output_dir, f"{model_name}.model")
compiled_model_path = os.path.join(output_dir, f"{model_name}_compiled.model")
⋮----
final_model_path = compiled_model_path
⋮----
final_model_path = model_path
⋮----
# 搜索checkpoints目录
checkpoints_dir = os.path.join(output_dir, "checkpoints")
⋮----
model_files = list(Path(checkpoints_dir).glob("*.model"))
⋮----
final_model_path = str(model_files[0])
⋮----
final_loss=None,  # 可以从日志中解析
⋮----
"""加载MACE计算器"""
⋮----
calculators = []
⋮----
calc = MACECalculator(
⋮----
"""创建平均力计算器 (用于ensemble)"""
⋮----
class AverageForceCalculator(Calculator)
⋮----
"""计算多个模型的平均力、能量和费米能级"""
implemented_properties = ['forces', 'energy', 'potential']
⋮----
def __init__(self, calculators, **kwargs)
⋮----
total_forces = 0
total_energy = 0
total_mu = 0
all_forces = []
all_mu = []
⋮----
atoms_copy = copy.deepcopy(atoms)
⋮----
n = len(self.calculators)
⋮----
# 计算不确定性
all_forces_array = np.array(all_forces)
force_std = np.std(all_forces_array, axis=0)
total_std = np.sqrt(np.sum(force_std**2, axis=1))
⋮----
def get_max_std(self)
⋮----
def get_mu_std(self)
⋮----
def get_mu(self, atoms=None)
⋮----
"""
        运行常电位分子动力学模拟

        Parameters
        ----------
        structure : str or ase.Atoms or pymatgen.Structure
            初始结构（带有吸附物的表面）。需要包含'electron'属性。
        model_paths : list of str
            CP-MACE模型文件路径列表。建议使用2个以上模型进行ensemble计算。
        target_potential : float
            目标电极电位 (V vs. SHE)
        initial_electron_number : float, optional
            初始电子数。如果不指定，从结构的info中读取。
        temperature : float, default=300.0
            模拟温度 (K)
        timestep : float, default=1.0
            时间步长 (fs)
        steps : int, default=1000
            模拟步数
        integrator : str, default="NoseHoover"
            积分器类型: "NoseHoover", "Langevin", "VelocityVerlet"
        ttime : float, default=40.0
            热浴耦合时间 (fs)
        eta_length : int, default=2
            Nose-Hoover链长度
        Mne : float, optional
            电子质量参数。如果不指定，使用初始电子数。
        constraints : list, optional
            约束列表，格式: [[type, atom1, atom2, distance], ...]
            type=0表示固定距离约束
        constraint_increment : float, default=0.0
            每步约束距离增量 (Å)，用于慢增长模拟
        read_velocity : bool, default=False
            是否从结构文件读取初始速度
        force_threshold : float, default=0.15
            力标准差阈值，超过此值收集结构 (eV/Å)
        fermi_threshold : float, default=0.04
            费米能级标准差阈值，超过此值收集结构 (V)
        save_frequency : int, default=1
            保存频率 (步数)
        output_dir : str, optional
            输出目录
        run_name : str, optional
            模拟名称

        Returns
        -------
        SimulationResult
            模拟结果
        """
⋮----
# 转换结构
atoms = self._convert_structure(structure)
⋮----
formula = atoms.get_chemical_formula()
⋮----
output_dir = os.path.join(
⋮----
# 设置电子数
⋮----
Mne = atoms.info['electron']
⋮----
# 设置dtype
⋮----
# 加载计算器
⋮----
calculators = self._load_calculators(model_paths)
⋮----
avg_calculator = self._create_average_calculator(calculators)
⋮----
# 创建积分器配置
integrator_config = {
⋮----
# 初始化速度
⋮----
# 创建积分器
⋮----
# 动态导入积分器
slow_growth_path = os.path.join(self.cp_mace_root, "simulation", "slow_growth")
⋮----
dyn = getattr(md_integrator, integrator)(atoms, **integrator_config)
⋮----
# 设置轨迹和日志
traj_path = os.path.join(output_dir, "atoms.traj")
traj = Trajectory(traj_path, 'w', atoms)
⋮----
# 运行模拟
⋮----
energy_history = []
fermi_history = []
electron_history = []
force_std_history = []
fermi_std_history = []
high_uncertainty_structures = []
trajectory = []
⋮----
iterator = tqdm(range(steps))
⋮----
iterator = range(steps)
⋮----
completed_steps = 0
⋮----
# 收集数据
energy = atoms.get_potential_energy()
⋮----
calc = atoms.calc
⋮----
fermi = calc.get_mu()
⋮----
fermi = calc.results['potential']
⋮----
fermi = atoms.info.get('potential', 0)
⋮----
# 不确定性
⋮----
force_std = calc.get_max_std()
fermi_std = calc.get_mu_std()
⋮----
force_std = 0
fermi_std = 0
⋮----
# 检查阈值，收集高不确定性结构
⋮----
# 清理GPU内存
⋮----
# 保存最终结构
final_path = os.path.join(output_dir, "final.xyz")
⋮----
# 保存高不确定性结构
⋮----
uncertain_path = os.path.join(output_dir, "high_uncertainty.xyz")
⋮----
result = SimulationResult(
⋮----
"""
        运行慢增长模拟（用于自由能计算）

        Parameters
        ----------
        structure : str or ase.Atoms
            初始结构
        model_paths : list of str
            模型路径列表
        target_potential : float
            目标电极电位 (V)
        atom1_index : int
            约束原子1的索引
        atom2_index : int
            约束原子2的索引
        initial_distance : float
            初始距离 (Å)
        final_distance : float
            最终距离 (Å)
        steps : int, default=1000
            模拟步数
        temperature : float, default=300.0
            温度 (K)
        timestep : float, default=1.0
            时间步长 (fs)
        output_dir : str, optional
            输出目录
        **kwargs
            其他参数传递给simulate()

        Returns
        -------
        SimulationResult
            模拟结果
        """
increment = (final_distance - initial_distance) / steps
⋮----
constraints = [[0, atom1_index, atom2_index, initial_distance]]
⋮----
def __del__(self)
⋮----
"""清理临时目录"""
⋮----
# 便捷函数
⋮----
"""
    便捷函数：训练CP-MACE模型

    Parameters
    ----------
    train_file : str
        训练数据文件路径
    cp_mace_root : str
        CP-MACE项目根目录
    model_name : str, default="CP_MACE_model"
        模型名称
    use_gpu : bool, default=True
        是否使用GPU
    output_dir : str, optional
        输出目录
    **kwargs
        其他参数传递给CPMACEPredictor.train()

    Returns
    -------
    TrainingResult
        训练结果
    """
predictor = CPMACEPredictor(
⋮----
"""
    便捷函数：运行常电位分子动力学模拟

    Parameters
    ----------
    structure : str or ase.Atoms
        初始结构
    model_paths : list of str
        模型路径列表
    cp_mace_root : str
        CP-MACE项目根目录
    target_potential : float
        目标电极电位 (V)
    temperature : float, default=300.0
        温度 (K)
    steps : int, default=1000
        模拟步数
    use_gpu : bool, default=True
        是否使用GPU
    output_dir : str, optional
        输出目录
    **kwargs
        其他参数

    Returns
    -------
    SimulationResult
        模拟结果
    """
⋮----
parser = argparse.ArgumentParser(
subparsers = parser.add_subparsers(dest="command", help="Commands")
⋮----
# Train command
train_parser = subparsers.add_parser("train", help="Train CP-MACE model")
⋮----
# Simulate command
sim_parser = subparsers.add_parser("simulate", help="Run CP-MD simulation")
⋮----
args = parser.parse_args()
⋮----
result = predictor.train(
⋮----
result = predictor.simulate(
````

## File: core/reconstruction/electrochemical_surface_predictor.py
````python
"""
Electrochemical Surface Predictor - CP-MACE + VSSR-MC 集成

这个模块实现了策略B：On-the-Fly完全集成
将CP-MACE作为VSSR-MC的计算器，实现真正的电化学模拟。

核心功能：
1. 使用CP-MACE在恒电势下计算能量（包含电场效应）
2. VSSR-MC进行表面结构采样（吸附/解吸/交换）
3. 基于Grand Potential的Metropolis判据（考虑pH和电势）
4. 支持Active Learning（可选）

使用方式:
    from electrochemical_surface_predictor import ElectrochemicalSurfacePredictor

    # 初始化预测器
    predictor = ElectrochemicalSurfacePredictor(
        cp_mace_root="/path/to/CP-MACE",
        surface_sampling_root="/path/to/surface-sampling",
        use_gpu=True,
    )

    # 在指定pH和电势下采样表面结构
    result = predictor.sample_surface(
        surface="slab.vasp",
        adsorbates=["O", "OH"],
        potential_she=1.0,  # V vs. SHE
        ph=12.0,
        temperature=300.0,
        total_sweeps=100,
    )

    print(result.summary())

作者: Claude
"""
⋮----
# 物理常数
KB = 8.617333262e-5  # eV/K
KB_J = 1.380649e-23  # J/K
E_CHARGE = 1.602176634e-19  # C
FARADAY = 96485.33212  # C/mol
R_GAS = 8.314462618  # J/(mol·K)
⋮----
# SHE到真空能级的转换 (常用值)
SHE_TO_VACUUM = -4.44  # eV (可根据CP-MACE训练基准调整)
⋮----
@dataclass
class ElectrochemicalConditions
⋮----
"""电化学条件"""
potential_she: float  # 电极电势 vs. SHE (V)
ph: float  # pH值
temperature: float  # 温度 (K)
⋮----
@property
    def potential_vacuum(self) -> float
⋮----
"""转换为真空能级 (eV)"""
⋮----
@property
    def kT(self) -> float
⋮----
"""热能 kT (eV)"""
⋮----
def __repr__(self)
⋮----
@dataclass
class ChemicalPotentials
⋮----
"""化学势数据"""
species: Dict[str, float]  # 各物种的化学势 (eV)
reference_state: str  # 参考态描述
⋮----
def get(self, species: str, default: float = 0.0) -> float
⋮----
@dataclass
class SampledConfiguration
⋮----
"""采样的表面配置"""
sweep_number: int
atoms: "ase.Atoms"
electronic_grand_potential: float  # Ω_el from CP-MACE (eV)
chemical_potential_correction: float  # Δμ_ions (eV)
total_grand_potential: float  # Ω_total (eV)
fermi_level: float  # 费米能级 (V)
electron_number: float  # 电子数
composition: Dict[str, int]  # 表面组成
num_adsorbates: int
acceptance_rate: float
force_uncertainty: float = 0.0  # 力的不确定性
fermi_uncertainty: float = 0.0  # 费米能级不确定性
⋮----
formula = self.atoms.get_chemical_formula() if self.atoms else "N/A"
⋮----
@dataclass
class ElectrochemicalSamplingResult
⋮----
"""电化学表面采样结果"""
surface_name: str
conditions: ElectrochemicalConditions
total_sweeps: int
sweep_size: int
configurations: List[SampledConfiguration]
grand_potential_history: List[float]
fermi_level_history: List[float]
electron_history: List[float]
acceptance_history: List[float]
best_configuration: SampledConfiguration
lowest_energy_configuration: SampledConfiguration
high_uncertainty_configs: List[SampledConfiguration] = field(default_factory=list)
output_dir: Optional[str] = None
⋮----
def get_top_n(self, n: int) -> List[SampledConfiguration]
⋮----
"""获取能量最低的前N个配置"""
sorted_configs = sorted(self.configurations, key=lambda x: x.total_grand_potential)
⋮----
def summary(self, top_n: int = 5) -> str
⋮----
"""生成结果摘要"""
lines = [
⋮----
class ElectrochemicalEvaluator
⋮----
"""
    电化学能量评估器

    封装CP-MACE计算器，并实现Grand Potential的计算。

    Grand Potential计算公式:
    Ω_total = Ω_el(U) + Δμ_ions(pH) + E_cap

    其中:
    - Ω_el(U): CP-MACE在恒电势下输出的电子巨势
    - Δμ_ions(pH): 离子化学势修正 (基于Nernst方程)
    - E_cap: 电容修正项 (可选)
    """
⋮----
# 标准化学势参考值 (eV)
# 这些值需要根据实际DFT计算校准
REFERENCE_ENERGIES = {
⋮----
"H2": -6.77,  # H2气体分子
"H2O": -14.22,  # H2O液态
"O2": -9.86,  # O2气体分子
⋮----
"""
        初始化电化学评估器

        Parameters
        ----------
        cp_mace_root : str
            CP-MACE根目录
        model_paths : list of str
            CP-MACE模型路径列表
        conditions : ElectrochemicalConditions
            电化学条件
        device : str
            计算设备
        use_ensemble : bool
            是否使用ensemble计算
        verbose : bool
            是否输出详细信息
        """
⋮----
# 添加路径
⋮----
# 延迟加载计算器
⋮----
def _load_calculators(self)
⋮----
"""加载CP-MACE计算器"""
⋮----
calc = MACECalculator(
⋮----
def _create_average_calculator(self)
⋮----
"""创建平均力计算器"""
⋮----
calculators = self._load_calculators()
⋮----
class CPMACEAverageCalculator(Calculator)
⋮----
"""CP-MACE Ensemble计算器"""
implemented_properties = ['energy', 'forces', 'potential']
⋮----
def __init__(self, calculators, **kwargs)
⋮----
total_forces = 0
total_energy = 0
total_potential = 0
all_forces = []
all_potentials = []
⋮----
atoms_copy = copy.deepcopy(atoms)
⋮----
n = len(self.calculators)
⋮----
# 计算不确定性
all_forces_array = np.array(all_forces)
force_std = np.std(all_forces_array, axis=0)
total_std = np.sqrt(np.sum(force_std**2, axis=1))
⋮----
def get_fermi_level(self)
⋮----
def get_force_uncertainty(self)
⋮----
def get_fermi_uncertainty(self)
⋮----
def get_calculator(self)
⋮----
"""获取计算器"""
⋮----
calcs = self._load_calculators()
⋮----
"""
        计算化学势修正项 Δμ_ions(pH)

        基于计算氢电极(CHE)模型:
        - 吸附O*:  ΔG = G(H2O) - G(H2) - 1/2*G(O2) - kT*ln(10)*pH
        - 吸附OH*: ΔG = G(H2O) - 1/2*G(H2) - kT*ln(10)*pH
        - 吸附H*:  ΔG = 1/2*G(H2) - e*U_SHE

        Parameters
        ----------
        added_species : dict
            添加/移除的物种计数, 如 {"O": 1, "H": -1}
            正数表示添加，负数表示移除

        Returns
        -------
        float
            化学势修正 (eV)
        """
correction = 0.0
kT = self.conditions.kT
pH = self.conditions.ph
U = self.conditions.potential_she
⋮----
# O* 吸附: 来自H2O，释放H2
# μ_O = μ_H2O - μ_H2 + kT*ln(10)*pH - eU
# 在pH=0, U=0时，μ_O ≈ 2.46 eV (相对于1/2 O2)
mu_O = (self.REFERENCE_ENERGIES["H2O"]
⋮----
# OH* 吸附: 来自H2O，释放1/2 H2
# μ_OH = μ_H2O - 0.5*μ_H2 + kT*ln(10)*pH - 0.5*eU
mu_OH = (self.REFERENCE_ENERGIES["H2O"]
⋮----
# H* 吸附: 来自溶液中的H+
# μ_H = 0.5*μ_H2 - kT*ln(10)*pH + eU
mu_H = (0.5 * self.REFERENCE_ENERGIES["H2"]
⋮----
# 其他元素: 使用默认化学势 (0)
# 用户可以通过继承此类来添加更多物种
⋮----
"""
        计算总巨势 Ω_total

        Ω_total = Ω_el(U) + Δμ_ions(pH)

        Parameters
        ----------
        atoms : ase.Atoms
            原子结构（需要包含'electron'属性）
        added_species : dict, optional
            相对于参考态添加/移除的物种
        relax : bool
            是否进行结构弛豫
        fmax : float
            弛豫收敛标准
        steps : int
            最大弛豫步数

        Returns
        -------
        tuple
            (Ω_total, Ω_el, Δμ_ions, fermi_level, electron_number)
        """
calc = self.get_calculator()
atoms = atoms.copy()
⋮----
# 确保有电子数属性
⋮----
# 估计初始电子数（可以基于元素的原子序数）
total_electrons = sum(atoms.get_atomic_numbers())
⋮----
# 可选：结构弛豫
⋮----
opt = BFGS(atoms, logfile=None)
⋮----
# 计算CP-MACE能量（电子巨势）
omega_el = atoms.get_potential_energy()
⋮----
# 获取费米能级
⋮----
fermi_level = calc.get_fermi_level()
⋮----
fermi_level = calc.results['potential']
⋮----
fermi_level = atoms.info.get('potential', 0)
⋮----
electron_number = atoms.info.get('electron', 0)
⋮----
# 计算化学势修正
⋮----
delta_mu = self.calculate_chemical_potential_correction(added_species)
⋮----
delta_mu = 0.0
⋮----
# 总巨势
omega_total = omega_el + delta_mu
⋮----
"""
        Metropolis判据

        P_accept = min(1, exp(-ΔΩ/kT))

        Parameters
        ----------
        delta_omega : float
            巨势变化 (eV)
        temperature_kT : float, optional
            温度 (kT单位)，默认使用conditions中的值

        Returns
        -------
        bool
            是否接受
        """
⋮----
kT = temperature_kT
⋮----
probability = np.exp(-delta_omega / kT)
⋮----
class ElectrochemicalSurfacePredictor
⋮----
"""
    电化学表面预测器 (CP-MACE + VSSR-MC)

    实现策略B：On-the-Fly完全集成
    将CP-MACE作为VSSR-MC的计算器，在恒电势下进行表面结构采样。

    Parameters
    ----------
    cp_mace_root : str
        CP-MACE项目根目录
    surface_sampling_root : str
        surface-sampling项目根目录
    model_paths : list of str
        CP-MACE模型路径列表
    device : str, default="cuda"
        计算设备
    cutoff : float, default=6.0
        截断半径 (Å)
    work_dir : str, optional
        工作目录
    keep_files : bool, default=True
        是否保留输出文件
    verbose : bool, default=True
        是否输出详细信息
    """
⋮----
# 设置模型路径
⋮----
# 使用CP-MACE自带的示例模型
default_models = [
⋮----
# 设置工作目录
⋮----
# 验证环境
⋮----
def _validate_environment(self)
⋮----
"""验证环境"""
⋮----
mace_init = os.path.join(self.cp_mace_root, "mace", "__init__.py")
⋮----
mcmc_init = os.path.join(self.surface_sampling_root, "mcmc", "__init__.py")
⋮----
def _log(self, message: str)
⋮----
"""输出日志"""
⋮----
def _convert_surface(self, surface) -> "ase.Atoms"
⋮----
"""转换输入结构"""
⋮----
"""
        生成虚拟吸附位点

        基于VSSR-MC的方法，在表面上方生成ontop/bridge/hollow位点
        """
⋮----
# 准备AtomsBatch
device = "cpu"  # 仅用于生成位点
slab_batch = get_atoms_batch(
⋮----
# 创建SurfaceSystem (不需要计算器)
system_settings = {
⋮----
# 使用dummy计算器
class DummyCalc
⋮----
def get_potential_energy(self, atoms=None)
def set(self, **kwargs)
⋮----
surface_system = SurfaceSystem(
⋮----
"""
        在指定电化学条件下采样表面结构

        Parameters
        ----------
        surface : str or ase.Atoms or pymatgen.Structure
            表面结构
        adsorbates : list of str
            可吸附的物种列表 (如 ["O", "OH", "H"])
        potential_she : float
            电极电势 vs. SHE (V)
        ph : float
            溶液pH值
        temperature : float, default=300.0
            温度 (K)
        initial_electron_number : float, optional
            初始电子数
        total_sweeps : int, default=100
            MC采样轮数
        sweep_size : int, default=20
            每轮采样步数
        canonical : bool, default=False
            是否使用canonical (NVT) 采样
        num_adsorbates : int, default=0
            固定吸附物数量 (canonical模式)
        force_threshold : float, default=0.15
            力不确定性阈值，用于Active Learning
        fermi_threshold : float, default=0.04
            费米能级不确定性阈值
        relax_each_step : bool, default=False
            是否在每步后进行结构弛豫
        output_dir : str, optional
            输出目录
        run_name : str, optional
            运行名称

        Returns
        -------
        ElectrochemicalSamplingResult
            采样结果
        """
⋮----
# 设置电化学条件
conditions = ElectrochemicalConditions(
⋮----
# 转换结构
slab = self._convert_surface(surface)
⋮----
surface_formula = slab.get_chemical_formula()
⋮----
# 设置输出目录
⋮----
output_dir = os.path.join(
⋮----
# 设置初始电子数
⋮----
# 估计初始电子数（基于原子的价电子）
total_electrons = sum(slab.get_atomic_numbers())
⋮----
# 创建电化学评估器
⋮----
evaluator = ElectrochemicalEvaluator(
⋮----
# 生成虚拟吸附位点
⋮----
ads_sites = self._generate_virtual_sites(slab, adsorbates)
⋮----
# 初始化MC采样
⋮----
# 计算初始巨势
current_atoms = slab.copy()
⋮----
# 初始化历史记录
configurations = []
grand_potential_history = []
fermi_history = []
electron_history = []
acceptance_history = []
high_uncertainty_configs = []
⋮----
# 追踪当前表面组成
current_composition = Counter(current_atoms.get_chemical_symbols())
adsorbate_counts = {ads: 0 for ads in adsorbates}
⋮----
# MC采样循环
total_accepted = 0
total_trials = 0
⋮----
sweep_iterator = tqdm(range(total_sweeps), desc="Sweeps")
⋮----
sweep_iterator = range(total_sweeps)
⋮----
sweep_accepted = 0
⋮----
# 选择试探步类型
# 1. 添加吸附物
# 2. 移除吸附物
# 3. 交换吸附物位置
move_type = np.random.choice(['add', 'remove', 'swap'],
⋮----
trial_atoms = current_atoms.copy()
added_species = {}
⋮----
# 添加吸附物
species = np.random.choice(adsorbates)
site_idx = np.random.randint(len(ads_sites))
site_pos = ads_sites[site_idx]
⋮----
# 添加原子
⋮----
# 移除吸附物（只能移除表面吸附的原子）
removable_indices = []
⋮----
# 检查是否在表面（z坐标较高）
⋮----
remove_idx = np.random.choice(removable_indices)
removed_species = trial_atoms[remove_idx].symbol
⋮----
# 没有可移除的原子，跳过
⋮----
# 交换吸附物位置
adsorbate_indices = []
⋮----
pos1 = trial_atoms[idx1].position.copy()
pos2 = trial_atoms[idx2].position.copy()
⋮----
# 保持电子数不变（从current_atoms继承）
⋮----
# 计算试探结构的巨势
⋮----
# Metropolis判据
delta_omega = trial_omega - current_omega
accepted = evaluator.metropolis_criterion(delta_omega, conditions.kT)
⋮----
current_atoms = trial_atoms
current_omega = trial_omega
current_omega_el = trial_omega_el
current_fermi = trial_fermi
current_ne = trial_ne
⋮----
# 更新组成
⋮----
# 记录这一轮的结果
acceptance_rate = sweep_accepted / sweep_size if sweep_size > 0 else 0
⋮----
# 获取不确定性
calc = evaluator.get_calculator()
force_uncertainty = 0.0
fermi_uncertainty = 0.0
⋮----
force_uncertainty = calc.get_force_uncertainty()
⋮----
fermi_uncertainty = calc.get_fermi_uncertainty()
⋮----
config = SampledConfiguration(
⋮----
# 检查不确定性，收集高不确定性结构 (Active Learning)
⋮----
# 找到最优结构
lowest_idx = np.argmin(grand_potential_history)
lowest_energy_config = configurations[lowest_idx]
best_config = configurations[-1]
⋮----
# 保存结果
⋮----
# 保存最终结构
⋮----
# 保存所有结构
all_atoms = [c.atoms for c in configurations]
⋮----
# 保存高不确定性结构
⋮----
uncertain_atoms = [c.atoms for c in high_uncertainty_configs]
⋮----
# 保存设置和历史
settings = {
⋮----
# 创建结果对象
result = ElectrochemicalSamplingResult(
⋮----
def __del__(self)
⋮----
"""清理临时目录"""
⋮----
# 便捷函数
⋮----
"""
    便捷函数：在指定电化学条件下采样表面结构

    Parameters
    ----------
    surface : str or ase.Atoms
        表面结构
    adsorbates : list of str
        可吸附物种
    potential_she : float
        电极电势 vs. SHE (V)
    ph : float
        溶液pH值
    cp_mace_root : str
        CP-MACE根目录
    surface_sampling_root : str
        surface-sampling根目录
    model_paths : list of str, optional
        模型路径
    temperature : float, default=300.0
        温度 (K)
    total_sweeps : int, default=100
        采样轮数
    device : str, default="cuda"
        计算设备
    output_dir : str, optional
        输出目录
    **kwargs
        其他参数

    Returns
    -------
    ElectrochemicalSamplingResult
        采样结果
    """
predictor = ElectrochemicalSurfacePredictor(
⋮----
parser = argparse.ArgumentParser(
⋮----
args = parser.parse_args()
⋮----
result = predictor.sample_surface(
````

## File: core/reconstruction/run_vssr_mc_subprocess.py
````python
#!/usr/bin/env python3
"""
Subprocess script for running VSSR-MC in isolation.

This script runs VSSR-MC in a separate process to avoid GPU/MPI state corruption
from other models (SurFF, AdsorbDiff) that may have been loaded in the main process.
"""
⋮----
# Increase recursion limit at module level
⋮----
# Add paths
_this_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.join(_this_dir, "..", "..")
sys.path.insert(0, _project_root)  # project root for `core.reconstruction.*`
sys.path.insert(0, os.path.join(_this_dir, ".."))  # core/
sys.path.insert(0, _this_dir)  # core/reconstruction/
⋮----
def patch_chgnet_cutoff()
⋮----
"""Monkey-patch CHGNet's atom_graph_cutoff to 10.0 to avoid isolated atom issues"""
⋮----
# Store original function
original_convert_data_batch = chgnet_io.convert_data_batch
⋮----
def patched_convert_data_batch(data_batch, cutoff=5.0, shuffle=True)
⋮----
"""Patched version that uses higher atom_graph_cutoff"""
⋮----
detached_batch = batch_detach(data_batch)
nxyz = detached_batch["nxyz"]
atoms_batch = AtomsBatch(
atoms_list = atoms_batch.get_list_atoms()
⋮----
pymatgen_structures = [AseAtomsAdaptor.get_structure(ab) for ab in atoms_list]
⋮----
energies = torch.atleast_1d(data_batch.get("energy"))
⋮----
energies_per_atom = energies
⋮----
energies_per_atom = torch.Tensor([0.0] * len(pymatgen_structures))
⋮----
energy_grads = data_batch.get("energy_grad")
num_atoms = detach(data_batch["num_atoms"]).tolist()
stresses = data_batch.get("stress")
magmoms = data_batch.get("magmoms")
⋮----
forces = [-x for x in energy_grads] if isinstance(energy_grads, list) else -energy_grads
⋮----
forces = None
⋮----
forces = torch.split(torch.atleast_2d(forces), num_atoms)
⋮----
forces = [torch.zeros_like(torch.Tensor(ab.get_positions())) for ab in atoms_list]
⋮----
stresses = torch.split(torch.atleast_2d(stresses), num_atoms)
⋮----
magmoms = torch.split(torch.atleast_2d(magmoms), num_atoms)
⋮----
# Use higher atom_graph_cutoff to handle edge cases
graph_converter = CrystalGraphConverter(
⋮----
atom_graph_cutoff=10.0,  # Increased from default 6.0
bond_graph_cutoff=5.0,   # Increased from default 3.0
⋮----
graph_converter=graph_converter,  # Pass custom graph_converter
⋮----
# Apply patch
⋮----
def main()
⋮----
parser = argparse.ArgumentParser(description="Run VSSR-MC in subprocess")
⋮----
args = parser.parse_args()
⋮----
# Load surface (with adsorbate)
⋮----
surface = pickle.load(f)
⋮----
# Load config
⋮----
config = json.load(f)
⋮----
# Load clean_slab (without adsorbate) for virtual site generation
clean_slab = None
⋮----
clean_slab = pickle.load(f)
⋮----
# Import VSSR-MC predictor
⋮----
# Create predictor
predictor = VSSRMCPredictor(
⋮----
# Run sampling
sample_kwargs = dict(
# Allow passing pre-computed offset_data to ensure consistent
# energy reference across multiple reconstruction steps
⋮----
result = predictor.sample(**sample_kwargs)
⋮----
# Save result - include all VSSRMCResult and SampledStructure fields
output_data = {
⋮----
error_data = {
````

## File: core/reconstruction/uma_surface_calculator.py
````python
"""
UMA Surface Calculator Adapters for VSSR-MC.

Wraps FairChem's FAIRChemCalculator (UMA model) to work with VSSR-MC's
SurfaceSystem, which requires a `surface_energy` property.

Two modes:
- UMASurfaceCalculator: thermal (gas-solid), same formula as EnsembleNFFSurface
- UMAPourbaixCalculator: electrochemical (liquid-solid), same formula as NFFPourbaix
"""
⋮----
logger = logging.getLogger(__name__)
⋮----
KB_EV = 8.617333262e-5  # eV/K
⋮----
class UMASurfaceCalculator(Calculator)
⋮----
"""Wraps any ASE Calculator (e.g. UMA FAIRChemCalculator) with surface energy.

    Surface energy formula (same as EnsembleNFFSurface, Du et al. 2023):
        E_surf = E_total - N_ref × E_bulk_ref
                 - Σ_i (excess_i × E_bulk_i)
                 - Σ_i (excess_i × μ_i)
    """
⋮----
implemented_properties = ("energy", "forces", "surface_energy")
⋮----
def set(self, **kwargs) -> dict
⋮----
changed = {}
⋮----
# Store remaining kwargs in parameters dict for compatibility
⋮----
atoms = self.atoms
⋮----
# Delegate energy/forces to base calculator
base_props = [p for p in properties if p in ("energy", "forces")]
⋮----
def get_potential_energy(self, atoms=None, **kwargs)
⋮----
def _compute_surface_energy(self, atoms: ase.Atoms) -> float
⋮----
"""Same formula as EnsembleNFFSurface.get_surface_energy()."""
e_total = self.get_potential_energy(atoms)
⋮----
ads_count = Counter(atoms.get_chemical_symbols())
bulk_energies = self.offset_data.get("bulk_energies", {})
stoics = self.offset_data.get("stoics", {})
ref_formula = self.offset_data.get("ref_formula", "")
ref_element = self.offset_data.get("ref_element", "")
⋮----
# Subtract bulk reference energy
bulk_ref_en = ads_count[ref_element] * bulk_energies.get(ref_formula, 0.0)
⋮----
s_ele = stoics.get(ele, 0)
s_ref = stoics.get(ref_element, 1)
excess = ads_count[ele] - s_ele / s_ref * ads_count[ref_element]
⋮----
surface_energy = e_total - bulk_ref_en
⋮----
# Subtract chemical potential deviation
⋮----
class UMAPourbaixCalculator(UMASurfaceCalculator)
⋮----
"""Wraps UMA with Pourbaix grand potential for electrochemical MC.

    Exact same formulas as NFFPourbaix (mcmc/calculators/calculators.py),
    including adsorbate corrections (OH ZPE-TS).

    Accepts PourbaixAtom objects (from mcmc.pourbaix.atoms) or dicts with
    the same keys (atom_std_state_energy, num_e, num_H, delta_G2_std, species_conc).
    """
⋮----
self.temp = KB_EV * temperature_k  # kT in eV
⋮----
changed = super().set(**kwargs)
⋮----
@staticmethod
    def _pa_attr(pa, key, default=0.0)
⋮----
"""Get attribute from PourbaixAtom object or dict."""
⋮----
"""Pourbaix grand potential = -(ΔG₁ + ΔG₂)."""
⋮----
def _get_delta_g1(self, atoms: ase.Atoms) -> float
⋮----
"""Dissociation free energy (same as NFFPourbaix.get_delta_G1).

        ΔG₁ = Σ(count × μ_atom_std_state) - E_slab + adsorbate_corrections
        """
e_slab = self.get_potential_energy(atoms)
⋮----
# Sum standard-state chemical potentials
atoms_count = Counter(atoms.get_chemical_symbols())
sum_chem_pots = 0.0
⋮----
# Adsorbate corrections (OH ZPE-TS), same logic as NFFPourbaix
⋮----
formula = Formula(atoms.get_chemical_formula())
⋮----
HO_diff = max(formula["H"] - formula["O"], 0)
⋮----
formula_dict_to_subtract = (F("H2O") * HO_diff).count()
formula_dict = formula.count()
formula_dict = {
formula = F.from_dict(formula_dict)
⋮----
def _get_delta_g2(self, atoms: ase.Atoms) -> float
⋮----
"""Electrochemical free energy (same as NFFPourbaix.get_delta_G2).

        ΔG₂ = Σ_atoms(ΔG₂_std - n_e·U - 2.3·n_H·kT·pH + kT·ln(a))
        """
delta_g2 = 0.0
⋮----
pa = self.pourbaix_atoms[symbol]
n_e = self._pa_attr(pa, "num_e", 0)
n_h = self._pa_attr(pa, "num_H", 0)
dg2_std = self._pa_attr(pa, "delta_G2_std", 0.0)
conc = self._pa_attr(pa, "species_conc", 1e-6)
delta_g2_non_std = (
````

## File: core/reconstruction/vssr_mc_predictor.py
````python
"""
VSSR-MC (Virtual Surface Site Relaxation - Monte Carlo) Predictor

This module provides a wrapper class for the surface-sampling (VSSR-MC) algorithm
to sample surface reconstructions across compositional and configurational spaces.

The VSSR-MC algorithm generates virtual adsorption sites above a pristine surface
and uses Monte Carlo sampling to explore different surface reconstructions by
adding/removing atoms at these sites.

Usage:
    from vssr_mc_predictor import VSSRMCPredictor

    # Initialize predictor
    predictor = VSSRMCPredictor(
        surface_sampling_root="/path/to/surface-sampling",
        model_type="CHGNetNFF",
    )

    # Run sampling - offset_data is automatically generated
    result = predictor.sample(
        surface=slab,  # ASE Atoms or path to structure file
        adsorbates=["O", "Sr", "Ti"],  # elements for surface reconstruction
        total_sweeps=100,
        temperature=1.0,
    )

    print(result.summary())
"""
⋮----
@dataclass
class SampledStructure
⋮----
"""Information about a sampled structure"""
sweep_number: int
atoms: "ase.Atoms"
energy: float  # Surface excess energy (eV) - relative to bulk reference
num_adsorbates: int
acceptance_rate: float
total_energy: Optional[float] = None  # CHGNet absolute total energy (eV)
⋮----
def __repr__(self)
⋮----
formula = self.atoms.get_chemical_formula() if self.atoms else "N/A"
⋮----
@dataclass
class VSSRMCResult
⋮----
"""Complete VSSR-MC sampling result"""
surface_name: str
model_type: str
total_sweeps: int
sweep_size: int
temperature: float
canonical: bool
structures: List[SampledStructure]
energy_history: List[float]
acceptance_history: List[float]
adsorption_count_history: List[int]
best_structure: SampledStructure
lowest_energy_structure: SampledStructure
offset_data: Optional[Dict] = None
output_dir: Optional[str] = None
⋮----
def get_structures_by_energy(self, top_n: int = 10) -> List[SampledStructure]
⋮----
"""Get structures sorted by energy (lowest first)"""
sorted_structures = sorted(self.structures, key=lambda x: x.energy)
⋮----
def get_unique_structures(self, energy_threshold: float = 0.01) -> List[SampledStructure]
⋮----
"""Get unique structures based on energy threshold"""
⋮----
unique = [sorted_structures[0]]
⋮----
def summary(self, top_n: int = 5) -> str
⋮----
"""Generate summary of sampling results"""
lines = [
⋮----
e_total_str = f", E_total={s.total_energy:.4f} eV" if s.total_energy is not None else ""
⋮----
class VSSRMCPredictor
⋮----
"""
    VSSR-MC Surface Reconstruction Sampler

    This class wraps the surface-sampling (VSSR-MC) algorithm to provide
    a simple interface for sampling surface reconstructions across
    compositional and configurational spaces.

    The predictor automatically generates offset_data from the input structure
    using bulk energy calculations from the NFF model.

    Parameters
    ----------
    surface_sampling_root : str
        Path to the surface-sampling repository root
    model_type : str, default="CHGNetNFF"
        Type of NFF model: "CHGNetNFF", "PaiNN", or "NffScaleMACE"
    model_paths : list of str, optional
        Paths to trained model checkpoints. If not provided and model_type
        is "CHGNetNFF", uses pre-trained CHGNet model.
    device : str, default="cuda"
        Device to use: "cuda" or "cpu"
    cutoff : float, optional
        Neighbor list cutoff. Uses model-specific default if not provided.
    work_dir : str, optional
        Working directory. Creates temp dir if not specified.
    keep_files : bool, default=True
        Whether to keep output files after sampling
    verbose : bool, default=True
        Whether to print progress information
    """
⋮----
DEFAULT_CUTOFFS = {
⋮----
# Electrochemical (Pourbaix) parameters
⋮----
# Set up work directory
⋮----
# Add surface-sampling to path
⋮----
# Validate environment
⋮----
# Lazy loading of models and calculators
⋮----
def _validate_environment(self)
⋮----
"""Validate the environment and dependencies"""
⋮----
mcmc_init = os.path.join(self.surface_sampling_root, "mcmc", "__init__.py")
⋮----
def _log(self, message: str)
⋮----
"""Log a message"""
⋮----
def _setup_logger(self, log_file: Path) -> logging.Logger
⋮----
"""Set up logger for MCMC"""
⋮----
def _load_models(self)
⋮----
"""Load NFF models"""
⋮----
device = get_final_device(self.device)
⋮----
# Both CHGNetNFF and NffScaleMACE support loading pre-trained models without paths
⋮----
# Use higher cutoff (8.0) to avoid isolated atom issues when adsorbates
# are slightly far from the surface (e.g., after diffusion)
model = load_model(
⋮----
cutoff=8.0,  # Increased from default 5.0 to handle edge cases
⋮----
"""Load surface energy calculator.

        Selects from 4 calculator types based on model_type and electrochemical params:
          CHGNet + thermal  → EnsembleNFFSurface
          CHGNet + Pourbaix → NFFPourbaix
          UMA    + thermal  → UMASurfaceCalculator
          UMA    + Pourbaix → UMAPourbaixCalculator
        """
is_pourbaix = (potential_she is not None and ph is not None)
⋮----
calc = self._load_uma_calculator(calc_settings, is_pourbaix, potential_she, ph)
⋮----
calc = self._load_nff_calculator(calc_settings, is_pourbaix, potential_she, ph)
⋮----
def _load_nff_calculator(self, calc_settings, is_pourbaix, potential_she, ph)
⋮----
"""Load NFF-based calculator (CHGNetNFF, NffScaleMACE, PaiNN)."""
⋮----
models = self._load_models()
⋮----
model_units = models[0].units if self.model_type != "PaiNN" else "kcal/mol"
⋮----
calc = NFFPourbaix(
pourbaix_settings = dict(calc_settings or {})
⋮----
# Auto-generate pourbaix_atoms if not provided:
# each element gets a default PourbaixAtom with bulk energy from offset_data
⋮----
# Use VSSR-MC's generate_pourbaix_atoms() to derive thermodynamically
# consistent Pourbaix parameters from pymatgen phase/Pourbaix diagrams.
⋮----
calc = EnsembleNFFSurface(
⋮----
def _load_uma_calculator(self, calc_settings, is_pourbaix, potential_she, ph)
⋮----
"""Load UMA (FairChem) based calculator.

        Uses the same import pattern as core/pathway/fairchem_predictor.py:
        local checkpoint → HuggingFace fallback.
        """
⋮----
# Import fairchem components (same pattern as fairchem_predictor.py)
⋮----
uma_model_name = "uma-s-1p1"
uma_model_path = self.model_paths[0] if self.model_paths else None
⋮----
# Try local checkpoints first
predictor = None
local_candidates = [
⋮----
predictor = load_predict_unit(
⋮----
# Fallback to HuggingFace
⋮----
predictor = pretrained_mlip.get_predict_unit(
⋮----
base_calc = FAIRChemCalculator(predictor, task_name="oc20")
⋮----
offset_data = (calc_settings or {}).get("offset_data")
chem_pots = (calc_settings or {}).get("chem_pots")
⋮----
pourbaix_atoms = (calc_settings or {}).get("pourbaix_atoms", {})
⋮----
# Use VSSR-MC's generate_pourbaix_atoms() for thermodynamically
# consistent parameters, then convert to dict format for UMA adapter
all_elems = set()
⋮----
# UMAPourbaixCalculator accepts both PourbaixAtom objects and dicts
pourbaix_atoms = self._generate_pourbaix_atoms(
temperature_k = 300.0
⋮----
temperature_k = calc_settings["temperature"] / 8.617333262e-5
calc = UMAPourbaixCalculator(
⋮----
calc = UMASurfaceCalculator(
⋮----
def _generate_pourbaix_atoms(self, elements, phi, pH)
⋮----
"""Generate PourbaixAtom objects using VSSR-MC's native Pourbaix module.

        Uses pymatgen PhaseDiagram + PourbaixDiagram to derive thermodynamically
        consistent parameters (num_e, num_H, delta_G2_std) for each element.
        """
⋮----
# Build phase diagram and Pourbaix diagram from pymatgen/MP data
⋮----
clean_elements = sorted(set(e for e in elements if len(e) <= 2 and e not in ("X",)))
⋮----
# Try Materials Project API
⋮----
mpr = MPRester()
⋮----
# Get phase diagram entries
entries = mpr.get_entries_in_chemsys(
pd = PhaseDiagram(entries)
⋮----
# Get Pourbaix entries
pbx_entries = mpr.get_pourbaix_entries(clean_elements)
pbx = PourbaixDiagram(pbx_entries)
⋮----
pourbaix_atoms = generate_pourbaix_atoms(pd, pbx, phi, pH, clean_elements)
⋮----
# Fallback: construct minimal PourbaixAtom from element data
⋮----
pa = {}
⋮----
ox_states = Element(elem).common_oxidation_states
ne = max(s for s in ox_states if s > 0) if any(s > 0 for s in ox_states) else 2
⋮----
ne = 2
⋮----
# O and H with known electrochemistry
⋮----
def _convert_surface(self, surface) -> "ase.Atoms"
⋮----
"""Convert input surface to ASE Atoms"""
⋮----
# Try pymatgen Structure
⋮----
"""Prepare AtomsBatch from ASE Atoms"""
⋮----
def _compute_bulk_energy(self, formula: str, per_atom: bool = True) -> float
⋮----
"""
        Compute bulk energy for a given formula using the current model
        (NFF for CHGNet/PaiNN/MACE, FairChem for UMA).

        Parameters
        ----------
        formula : str
            Chemical formula (e.g., "Cu", "O2", "SrTiO3")
        per_atom : bool
            If True, return energy per atom

        Returns
        -------
        float
            Bulk energy in eV
        """
⋮----
atoms = None
⋮----
# Try to get bulk structure from Materials Project
⋮----
docs = mpr.materials.summary.search(
⋮----
docs_sorted = sorted(docs, key=lambda x: x.energy_per_atom if x.energy_per_atom else 0)
structure = docs_sorted[0].structure
atoms = AseAtomsAdaptor.get_atoms(structure)
⋮----
# Fall back to simple bulk structures for elements
⋮----
comp = Composition(formula)
⋮----
atoms = ase_bulk(formula)
⋮----
energy = self._compute_energy_uma(atoms)
⋮----
energy = self._compute_energy_nff(atoms)
⋮----
def _compute_energy_nff(self, atoms) -> float
⋮----
"""Compute energy using NFF model (CHGNet/PaiNN/MACE)."""
⋮----
slab_batch = self._prepare_slab_batch(atoms, device)
energy = calc.get_potential_energy(atoms=slab_batch)
⋮----
energy = float(energy.item() if hasattr(energy, 'item') else energy[0])
⋮----
def _compute_energy_uma(self, atoms) -> float
⋮----
"""Compute energy using UMA (FairChem) model."""
⋮----
# Reuse cached predictor if available
⋮----
calc = FAIRChemCalculator(self._uma_predictor, task_name="oc20")
atoms_copy = atoms.copy()
⋮----
"""
        Automatically generate offset_data from the slab composition.

        The offset_data contains:
        - bulk_energies: bulk energy per formula unit for each element and reference compound
        - stoics: stoichiometry of each element in the reference formula
        - ref_formula: reference compound formula (derived from slab)
        - ref_element: reference element for normalization (usually the most common cation)

        Parameters
        ----------
        slab : ase.Atoms
            The pristine slab structure
        adsorbates : list of str, optional
            Additional adsorbate elements

        Returns
        -------
        dict
            offset_data dictionary
        """
⋮----
# Get slab composition
slab_symbols = slab.get_chemical_symbols()
composition = Counter(slab_symbols)
elements = list(composition.keys())
⋮----
# Include adsorbate elements if provided
⋮----
# Determine reference formula from slab composition
# Use pymatgen to get reduced formula
comp = Composition(composition)
ref_formula = comp.reduced_formula
# Get stoichiometry from reduced formula, not original composition
reduced_comp = Composition(ref_formula)
el_amt = reduced_comp.get_el_amt_dict()
stoics = {str(el): int(amt) for el, amt in el_amt.items()}
⋮----
# Determine reference element (usually a cation that's not O/N/S/F/Cl)
# Priority: transition metals > main group metals > non-metals
anions = {"O", "N", "S", "F", "Cl", "Br", "I"}
cations = [el for el in stoics.keys() if el not in anions]
⋮----
ref_element = cations[0]  # First non-anion element
⋮----
ref_element = list(stoics.keys())[0]
⋮----
# Calculate bulk energies
bulk_energies = {}
⋮----
# First, try to compute bulk energy for the reference compound
⋮----
ref_energy = self._compute_bulk_energy(ref_formula, per_atom=False)
# Scale to per formula unit
n_atoms_per_fu = sum(stoics.values())
atoms_in_calc = len(slab)  # Approximate
⋮----
# Compute bulk energies for individual elements
⋮----
el_energy = self._compute_bulk_energy(el, per_atom=True)
⋮----
offset_data = {
⋮----
"""
        Generate default chemical potentials.

        By default, sets all chemical potentials to 0, which corresponds to
        elemental reference states.

        Parameters
        ----------
        slab : ase.Atoms
            The slab structure
        adsorbates : list of str, optional
            Adsorbate elements
        ref_element : str, optional
            Reference element (will have chem_pot = 0)

        Returns
        -------
        dict
            Chemical potentials for each element
        """
composition = Counter(slab.get_chemical_symbols())
⋮----
# Default: all chemical potentials = 0
chem_pots = {el: 0.0 for el in elements}
⋮----
"""
        Perform VSSR-MC sampling on a surface.

        This method samples surface reconstructions across compositional and
        configurational spaces using Monte Carlo with virtual adsorption sites.

        Parameters
        ----------
        surface : str or ase.Atoms or pymatgen.Structure
            Input surface structure (with adsorbate if doing canonical sampling)
        adsorbates : list of str, optional
            List of elements that can be added/removed at virtual sites
        canonical : bool, default=False
            If True, perform canonical (NVT) sampling with fixed composition.
            If False, perform semi-grand canonical sampling with variable composition.
        num_adsorbates : int, default=0
            Number of adsorbate atoms (required if canonical=True)
        total_sweeps : int, default=100
            Number of MC sweeps to perform
        sweep_size : int, default=20
            Number of MC steps per sweep
        temperature : float, default=1.0
            Temperature in kT units
        perform_annealing : bool, default=False
            Whether to perform simulated annealing
        annealing_alpha : float, default=0.99
            Annealing schedule parameter (T_new = alpha * T_old)
        chem_pots : dict, optional
            Chemical potentials for each element. If None, auto-generated with all 0.
        offset_data : dict, optional
            Offset data for surface energy calculation. If None and auto_offset=True,
            automatically generated from slab composition.
        auto_offset : bool, default=True
            If True and offset_data is None, automatically generate offset_data.
        calc_settings : dict, optional
            Additional calculator settings
        system_settings : dict, optional
            Additional system settings (surface_depth, near_reduce, etc.)
        output_dir : str, optional
            Directory to save results
        run_name : str, optional
            Name for this sampling run
        clean_slab : ase.Atoms, optional
            Clean slab without adsorbates, used for virtual site generation.
            If provided, virtual sites are generated on this clean surface
            (preventing sites from being placed above adsorbates).

        Returns
        -------
        VSSRMCResult
            Sampling results including all structures and statistics
        """
⋮----
# Convert and prepare surface
⋮----
surface_atoms = self._convert_surface(surface)
surface_formula = surface_atoms.get_chemical_formula()
⋮----
# Set up output directory
⋮----
output_dir = os.path.join(
⋮----
# Set up logger
logger = self._setup_logger(Path(output_dir) / "mc.log")
⋮----
# Auto-generate offset_data if needed
⋮----
offset_data = self._auto_generate_offset_data(surface_atoms, adsorbates)
⋮----
# Auto-generate chemical potentials if needed
⋮----
ref_element = offset_data.get("ref_element") if offset_data else None
chem_pots = self._auto_generate_chem_pots(surface_atoms, adsorbates, ref_element)
⋮----
# Prepare calculator settings
_calc_settings = calc_settings.copy() if calc_settings else {}
⋮----
# Prepare system settings
_system_settings = {
⋮----
"surface_depth": 2,  # Allow top 2 layers to move during reconstruction
⋮----
"planar_distance": 1.5,  # 虚拟位点在表面上方的高度
⋮----
# Inject external surface/adsorbate indices (overrides z-based detection)
⋮----
# 如果指定了max_site_height，使用它来覆盖planar_distance
# 这可以防止在吸附物上方生成虚拟位点
⋮----
max_height = _system_settings.pop("max_site_height")
⋮----
# Get device
⋮----
# Load calculator — select mode based on electrochemical params
eff_potential = potential_she if potential_she is not None else self.potential_she
eff_ph = ph if ph is not None else self.ph
⋮----
calculator = self._load_calculator(
⋮----
# 预先生成虚拟位点坐标（如果提供了clean_slab）
# 使用clean_slab（不含吸附物）生成虚拟位点，避免在吸附物上方生成
ads_coords = None
⋮----
pmg_clean = AseAtomsAdaptor.get_structure(clean_slab)
site_finder = AdsorbateSiteFinder(pmg_clean)
⋮----
ads_site_type = _system_settings.get("ads_site_type", "all")
ads_positions = site_finder.find_adsorption_sites(
⋮----
ads_coords = ads_positions
⋮----
# Prepare AtomsBatch (with adsorbate)
⋮----
slab_batch = self._prepare_slab_batch(surface_atoms, device)
⋮----
# Create SurfaceSystem
# 如果提供了ads_coords，SurfaceSystem将使用预先生成的虚拟位点
surface_system = SurfaceSystem(
⋮----
ads_coords=ads_coords,  # 使用在clean_slab上生成的虚拟位点
⋮----
# Save initial structure with virtual adsorption sites
⋮----
initial_energy = surface_system.get_surface_energy()
⋮----
# Set up MCMC
sampling_settings = {
⋮----
mcmc = MCMC(
⋮----
# Run sampling
⋮----
results = mcmc.run(
⋮----
# Process results
structures = []
⋮----
atoms = surf_sys.real_atoms.copy()
⋮----
atoms = surf_sys.atoms.copy()
⋮----
atoms = surf_sys
⋮----
# Compute absolute total energy from surface energy + bulk reference
total_e = None
⋮----
bulk_energies = offset_data.get("bulk_energies", {})
stoics = offset_data.get("stoics", {})
ref_element = offset_data.get("ref_element")
⋮----
atom_counts = Counter(atoms.get_chemical_symbols())
n_ref = atom_counts.get(ref_element, 0)
e_bulk_ref = bulk_energies[ref_element]
total_e = float(energy) + n_ref * e_bulk_ref
⋮----
sampled = SampledStructure(
⋮----
# Find best structures
lowest_energy_idx = np.argmin(results["energy_hist"])
lowest_energy_structure = structures[lowest_energy_idx]
⋮----
# Final structure
best_structure = structures[-1]
⋮----
# Save structures
⋮----
# Save offset_data
⋮----
# Save settings
all_settings = {
⋮----
result = VSSRMCResult(
⋮----
def __del__(self)
⋮----
"""Clean up temporary directory"""
⋮----
"""
    Convenience function for VSSR-MC surface reconstruction sampling.

    Parameters
    ----------
    surface : str or ase.Atoms
        Input surface structure
    surface_sampling_root : str
        Path to surface-sampling repository
    model_type : str, default="CHGNetNFF"
        Model type to use
    adsorbates : list of str, optional
        List of adsorbate elements
    total_sweeps : int, default=100
        Number of MC sweeps
    temperature : float, default=1.0
        Temperature in kT
    output_dir : str, optional
        Output directory
    device : str, default="cuda"
        Device to use
    **kwargs
        Additional arguments

    Returns
    -------
    VSSRMCResult
        Sampling results
    """
predictor = VSSRMCPredictor(
⋮----
parser = argparse.ArgumentParser(description="VSSR-MC Surface Reconstruction Sampler")
⋮----
args = parser.parse_args()
⋮----
result = predictor.sample(
````

## File: core/surface/__init__.py
````python

````

## File: core/surface/surff_predictor.py
````python
"""
SurFF Predictor - 金属间化合物表面能预测封装类

这个模块提供了一个完整的封装类，用于预测任意bulk结构的表面暴露情况。

使用方式:
    from surff_predictor import SurFFPredictor

    # 初始化预测器
    predictor = SurFFPredictor(
        surff_root="/path/to/SurFF",
        checkpoint_path="/path/to/checkpoint.pt",  # 可选，默认使用SurFF自带的模型
        use_gpu=True,
    )

    # 预测表面暴露
    results = predictor.predict(
        structure="POSCAR",  # 可以是文件路径、ASE Atoms或pymatgen Structure
        top_n=5,  # 返回最可能暴露的前N个表面
    )

作者: Claude
"""
⋮----
@dataclass
class SlabInfo
⋮----
"""表面slab信息"""
slab_id: str
miller_index: Tuple[int, int, int]
shift: float
num_atoms: int
formula: str
⋮----
@dataclass
class SurfaceResult
⋮----
"""单个表面的预测结果"""
⋮----
surface_energy: float  # eV/Å²
area_fraction: float   # 暴露面积比例 (0-1)
exposure_level: str    # "high", "medium", "low"
⋮----
def __repr__(self)
⋮----
@dataclass
class PredictionResult
⋮----
"""完整的预测结果"""
crystal_id: str
⋮----
total_slabs: int
surfaces: List[SurfaceResult]
wulff_shape_path: Optional[str] = None
raw_results: Dict = field(default_factory=dict)
⋮----
def get_top_n(self, n: int) -> List[SurfaceResult]
⋮----
"""获取暴露面积最大的前N个表面"""
sorted_surfaces = sorted(self.surfaces, key=lambda x: x.area_fraction, reverse=True)
⋮----
def get_exposed_surfaces(self, threshold: float = 0.01) -> List[SurfaceResult]
⋮----
"""获取所有暴露面积大于阈值的表面"""
⋮----
def to_dataframe(self) -> pd.DataFrame
⋮----
"""转换为DataFrame"""
data = []
⋮----
def summary(self, top_n: int = 5) -> str
⋮----
"""生成预测结果摘要"""
lines = [
⋮----
class SurFFPredictor
⋮----
"""
    SurFF表面能预测器

    这个类封装了SurFF项目的所有功能，提供简洁的API用于预测任意bulk结构的表面暴露情况。

    Parameters
    ----------
    surff_root : str
        SurFF项目的根目录路径
    checkpoint_path : str, optional
        模型checkpoint文件路径。如果不指定，将使用SurFF自带的模型
    use_gpu : bool, default=True
        是否使用GPU进行推理
    max_miller_index : int, default=2
        生成slab时的最大Miller指数
    min_slab_size : float, default=10.0
        最小slab厚度 (Å)
    min_vacuum_size : float, default=15.0
        最小真空层厚度 (Å)
    fix_distance : float, default=3.0
        固定原子的距离阈值 (Å)，距离表面超过此距离的原子将被固定
    relaxation_steps : int, default=200
        最大弛豫步数
    relaxation_fmax : float, default=0.05
        弛豫收敛力阈值 (eV/Å)
    max_step : float, default=0.03
        LBFGS优化器的最大步长
    work_dir : str, optional
        工作目录，用于存储中间文件。如果不指定，将创建临时目录
    keep_files : bool, default=False
        是否保留中间文件
    verbose : bool, default=True
        是否输出详细信息
    """
⋮----
# 设置checkpoint路径
⋮----
# 使用默认的checkpoint
default_ckpt = os.path.join(
⋮----
# 尝试ocp目录下的checkpoint
alt_ckpt = os.path.join(
⋮----
# 设置工作目录
⋮----
# 添加SurFF路径到系统路径
⋮----
# 验证环境
⋮----
def _validate_environment(self)
⋮----
"""验证运行环境"""
# 检查必要的目录
⋮----
# 检查必要的配置文件
config_path = os.path.join(self.ocp_root, "configs/equiformer_v2_002_relax.yml")
⋮----
def _log(self, message: str)
⋮----
"""输出日志"""
⋮----
def _convert_structure(self, structure) -> "pymatgen.core.Structure"
⋮----
"""将输入结构转换为pymatgen Structure"""
⋮----
# 文件路径
⋮----
# 尝试作为ASE Atoms处理
⋮----
"""生成所有可能的slab结构"""
⋮----
# 获取conventional cell
structure = SpacegroupAnalyzer(structure, symprec=0.01).get_conventional_standard_structure()
formula = structure.formula
⋮----
# 生成所有slabs
⋮----
slabs = generate_all_slabs(
⋮----
# 保存slab和bulk结构
slab_info = []
⋮----
miller_index = slab.miller_index
miller_str = ''.join(map(str, miller_index))
shift = slab.shift
num_atoms = len(slab.sites)
⋮----
slab_id = f"{crystal_id}_{j}"
⋮----
# 保存slab (带selective dynamics)
slab_path = os.path.join(slab_dir, slab_id)
poscar = self._fix_middle_atoms(slab)
⋮----
# 保存oriented unit cell
bulk_path = os.path.join(bulk_dir, slab_id)
ouc = slab.oriented_unit_cell
ouc_lattice = Lattice.from_parameters(
ouc = Structure(ouc_lattice, ouc.species, ouc.frac_coords, coords_are_cartesian=False)
⋮----
def _fix_middle_atoms(self, structure: "pymatgen.core.Structure") -> "pymatgen.io.vasp.Poscar"
⋮----
"""固定slab中间的原子"""
⋮----
selective_dynamics = []
⋮----
A = structure.lattice.matrix[0]
B = structure.lattice.matrix[1]
N = np.cross(A, B)
N = N / np.linalg.norm(N)
⋮----
distances = [abs(np.dot(N, site)) for site in structure.cart_coords]
d_max = max(distances)
d_min = min(distances)
⋮----
dist_from_top = d_max - d
dist_from_bottom = d - d_min
⋮----
def _generate_lmdb(self, slab_dir: str, lmdb_path: str, slab_info: pd.DataFrame)
⋮----
"""生成LMDB数据集"""
⋮----
slab_files = [f for f in os.listdir(slab_dir) if not f.endswith('.csv')]
dataset = [{'slab_id': f, 'POSCAR_pth': os.path.join(slab_dir, f)} for f in slab_files]
⋮----
db = lmdb.open(lmdb_path, map_size=10**9, subdir=False, meminit=False, map_async=True)
txn = db.begin(write=True)
⋮----
db_idx = 0
⋮----
poscar_path = item['POSCAR_pth']
sid = item['slab_id']
⋮----
crystal = Poscar.from_file(poscar_path).structure
natoms = len(crystal)
⋮----
pos = torch.Tensor(crystal.cart_coords.tolist())
⋮----
tags = crystal.site_properties['selective_dynamics']
tags = [1 if x[0] else 0 for x in tags]
tags = torch.LongTensor(tags)
⋮----
tags = torch.LongTensor([1] * natoms)
⋮----
fixed = (tags == 0).float()
lattice = torch.Tensor([crystal.lattice.matrix.tolist()])
atom_fea = torch.Tensor([crystal[i].specie.number for i in range(natoms)])
⋮----
data = Data(
⋮----
# 存储长度
⋮----
def _run_relaxation(self, lmdb_dir: str, traj_dir: str)
⋮----
"""运行MLFF弛豫"""
⋮----
# Convert absolute paths to relative paths from ocp root
lmdb_rel = os.path.relpath(os.path.abspath(lmdb_dir), self.ocp_root)
traj_rel = os.path.relpath(os.path.abspath(traj_dir), self.ocp_root)
⋮----
cmd = [
⋮----
result = subprocess.run(
⋮----
error_msg = result.stderr if result.stderr else "Unknown error"
⋮----
def _collect_results(self, traj_dir: str) -> pd.DataFrame
⋮----
"""收集弛豫结果"""
⋮----
traj_files = [f for f in os.listdir(traj_dir) if f.endswith('.traj')]
⋮----
results = []
⋮----
slab_id = os.path.splitext(traj_file)[0]
crystal_id = slab_id.split('_')[0]
⋮----
traj = Trajectory(os.path.join(traj_dir, traj_file))
final_atoms = traj[-1]
⋮----
energy = final_atoms.get_potential_energy()
cell = final_atoms.get_cell()
area = np.linalg.norm(np.cross(cell[0], cell[1]))
surface_energy = energy / (2 * area)
⋮----
"""计算Wulff构型"""
⋮----
# 合并结果
df = pd.merge(results_df, slab_info, on='slab_id', how='inner')
df = df.rename(columns={'surface_energy': 'surface_energy_pred'})
⋮----
# 读取晶体结构
crystal_path = os.path.join(crystal_dir, crystal_id)
crystal_struct = Poscar.from_file(crystal_path).structure
crystal_struct = SpacegroupAnalyzer(crystal_struct).get_conventional_standard_structure()
⋮----
# 解析Miller指数
def str_to_tuple(s)
⋮----
matches = re.findall(r'-?\d', str(s))
⋮----
miller_indices = [str_to_tuple(m) for m in df['miller_index'].values]
surface_energies = df['surface_energy_pred'].values.tolist()
⋮----
# 计算Wulff形状
wulff = WulffShape(crystal_struct.lattice, miller_indices, surface_energies)
⋮----
# 提取结果
wulff_results = {
⋮----
# 计算每个slab的暴露面积
⋮----
area = 0.0
⋮----
area = wulff.color_area[j] / sum(wulff.color_area) if sum(wulff.color_area) > 0 else 0
⋮----
# 保存Wulff形状图
wulff_shape_path = None
⋮----
matplotlib.use('Agg')  # 使用非交互式后端
⋮----
ax_3d = wulff.get_plot(
wulff_shape_path = os.path.join(save_dir, f"{crystal_id}_wulff.png")
⋮----
"""
        预测bulk结构的表面暴露情况

        Parameters
        ----------
        structure : str or pymatgen.Structure or ase.Atoms
            输入的bulk结构，可以是:
            - POSCAR/VASP格式文件路径
            - pymatgen Structure对象
            - ASE Atoms对象
        crystal_id : str, default="crystal"
            晶体ID，用于标识输出文件
        top_n : int, default=5
            返回最可能暴露的前N个表面
        save_wulff_shape : bool, default=True
            是否保存Wulff形状图
        output_dir : str, optional
            输出目录。如果不指定，使用工作目录

        Returns
        -------
        PredictionResult
            预测结果，包含所有表面的信息和Wulff形状
        """
# 设置输出目录
⋮----
output_dir = os.path.join(self.work_dir, crystal_id)
⋮----
# 创建子目录
crystal_dir = os.path.join(output_dir, "crystal")
slab_dir = os.path.join(output_dir, "slabs")
bulk_dir = os.path.join(output_dir, "bulks")
lmdb_dir = os.path.join(output_dir, "lmdb")
traj_dir = os.path.join(output_dir, "traj")
wulff_dir = os.path.join(output_dir, "wulff") if save_wulff_shape else None
⋮----
# 1. 转换并保存输入结构
⋮----
pmg_structure = self._convert_structure(structure)
formula = pmg_structure.formula
⋮----
# 保存晶体结构
⋮----
# 2. 生成slab结构
slab_info = self._generate_slabs(pmg_structure, crystal_id, slab_dir, bulk_dir)
⋮----
# 3. 生成LMDB数据集
lmdb_path = os.path.join(lmdb_dir, "relaxation.lmdb")
⋮----
# 4. 运行弛豫
⋮----
# 5. 收集结果
results_df = self._collect_results(traj_dir)
⋮----
# 6. 计算Wulff构型
⋮----
# 7. 构建返回结果
surfaces = []
⋮----
area = wulff_results['area_fractions'][i]
⋮----
level = "high"
⋮----
level = "medium"
⋮----
level = "low"
⋮----
result = PredictionResult(
⋮----
# 清理临时文件（如果需要）
⋮----
pass  # 保留输出目录的结果
⋮----
"""
        批量预测多个bulk结构的表面暴露情况

        Parameters
        ----------
        structures : list
            输入的bulk结构列表
        crystal_ids : list of str, optional
            晶体ID列表。如果不指定，将自动生成
        top_n : int, default=5
            返回最可能暴露的前N个表面
        save_wulff_shape : bool, default=True
            是否保存Wulff形状图
        output_dir : str, optional
            输出目录

        Returns
        -------
        list of PredictionResult
            每个结构的预测结果列表
        """
⋮----
crystal_ids = [f"crystal_{i}" for i in range(len(structures))]
⋮----
result = self.predict(
⋮----
def __del__(self)
⋮----
"""清理临时目录"""
⋮----
# 便捷函数
⋮----
"""
    便捷函数：预测bulk结构的表面暴露情况

    Parameters
    ----------
    structure : str or pymatgen.Structure or ase.Atoms
        输入的bulk结构
    surff_root : str
        SurFF项目根目录
    checkpoint_path : str, optional
        模型checkpoint路径
    top_n : int, default=5
        返回最可能暴露的前N个表面
    use_gpu : bool, default=True
        是否使用GPU
    output_dir : str, optional
        输出目录
    **kwargs
        其他参数传递给SurFFPredictor

    Returns
    -------
    PredictionResult
        预测结果
    """
predictor = SurFFPredictor(
⋮----
# 示例用法
⋮----
parser = argparse.ArgumentParser(description="SurFF Surface Exposure Predictor")
⋮----
args = parser.parse_args()
⋮----
result = predictor.predict(
⋮----
# 打印结果
⋮----
# 保存CSV
⋮----
csv_path = os.path.join(args.output_dir, f"{args.crystal_id}_results.csv")
````

## File: core/viz/__init__.py
````python

````

## File: core/viz/catalyst_surface_visualizer.py
````python
#!/usr/bin/env python
"""Catalyst Surface Visualizer with OVITO high-quality rendering only."""
⋮----
# Ensure Qt platform plugins are discoverable for pip-installed PySide6/OVITO.
_pyside6_plugins = os.path.join(
⋮----
# Optional imports for OVITO rendering (lazy-loaded to avoid startup instability).
OVITO_AVAILABLE = False
TACHYON_AVAILABLE = False
OSPRAY_AVAILABLE = False
OPENGL_AVAILABLE = False
ANARI_AVAILABLE = False
_OVITO_PROBED = False
⋮----
def _probe_ovito_runtime() -> None
⋮----
"""Probe OVITO runtime capabilities lazily and cache the result."""
⋮----
OVITO_AVAILABLE = True
⋮----
_OVITO_PROBED = True
⋮----
TACHYON_AVAILABLE = True
⋮----
OSPRAY_AVAILABLE = True
⋮----
OPENGL_AVAILABLE = True
⋮----
ANARI_AVAILABLE = True
⋮----
class CatalystSurfaceVisualizer
⋮----
"""High-quality visualizer for catalyst surfaces using OVITO renderers."""
⋮----
"""
        Initialize the visualizer.

        Args:
            elevation: Viewing angle elevation in degrees (default: 30)
            azimuth: Viewing angle azimuth in degrees (default: 45)
            scale: Scale factor for atom sizes (default: 1.0)
            quality: Rendering quality ('low', 'medium', 'high', 'ultra')
            background_color: Background color
            color_scheme: Color scheme ('jmol', 'cpk')
            renderer: Rendering backend ('tachyon', 'ospray', 'opengl', 'anari')
            show_cell: Show unit cell boundaries (default: False)
            auto_expand: Automatically expand small structures (default: True)
            expand_threshold: Minimum number of atoms before expansion (default: 15)
            supercell: Supercell dimensions for expansion (default: 2x2x1)
        """
⋮----
# Check renderer availability
ovito_renderers = {
⋮----
# Prefer isolated rendering via `ovitos` subprocess to avoid importing
# GUI/Jupyter extension stack in the main workflow process.
use_direct = os.environ.get("CATDT_OVITO_DIRECT", "0") == "1"
ovitos_exec = os.environ.get("CATDT_OVITOS_EXECUTABLE", "").strip() or shutil.which("ovitos")
⋮----
# Quality settings
quality_settings = {
⋮----
settings = quality_settings.get(quality, quality_settings['high'])
⋮----
def expand_structure_if_needed(self, atoms: Atoms) -> Tuple[Atoms, bool]
⋮----
"""
        Expand structure if it has fewer atoms than threshold.

        Args:
            atoms: ASE Atoms object

        Returns:
            (expanded_atoms, was_expanded): Tuple of atoms and whether expansion occurred
        """
⋮----
# Create supercell matrix
P = np.diag(self.supercell)
expanded = make_supercell(atoms, P)
⋮----
def get_atom_color(self, atomic_number: int) -> np.ndarray
⋮----
"""Get color for atomic number from jmol scheme."""
⋮----
def get_atom_radius(self, atomic_number: int) -> float
⋮----
"""Get radius for atomic number."""
⋮----
"""
        Render structure using OVITO renderers.

        By default, rendering is isolated in a subprocess to avoid hard
        process crashes (segfault) from OVITO/Tachyon C++ backends.
        """
⋮----
"""Isolate OVITO rendering in a subprocess for robustness."""
⋮----
interpreter = ovitos_exec if ovitos_exec else sys.executable
cmd = [
⋮----
env = os.environ.copy()
⋮----
timeout_s = int(str(os.environ.get("CATDT_OVITO_RENDER_TIMEOUT_SEC", "180")).strip() or "180")
⋮----
result = subprocess.run(
⋮----
err = (result.stderr or result.stdout or "").strip()
⋮----
err = err[-2000:]
⋮----
"""
        Render structure using OVITO renderers.

        Args:
            atoms: ASE Atoms object
            output_file: Output filename
            renderer_type: 'tachyon', 'ospray', 'opengl', or 'anari'
        """
⋮----
# Save atoms to temporary file for OVITO to import
⋮----
pipeline = None
data = None
viewport = None
renderer = None
⋮----
# Write atoms to temporary xyz file
⋮----
# Import with OVITO
pipeline = import_file(temp_file)
⋮----
# Add pipeline to scene for rendering
⋮----
# Compute pipeline
data = pipeline.compute()
⋮----
# Configure cell display
⋮----
# Control cell visualization based on user preference
⋮----
pass  # Cell visualization may not be available
⋮----
# Create viewport
viewport = Viewport()
⋮----
# Set camera position
positions = atoms.get_positions()
center = positions.mean(axis=0)
bbox_size = positions.max(axis=0) - positions.min(axis=0)
max_dim = np.max(bbox_size)
camera_distance = max_dim * 3.0
⋮----
# Camera position in spherical coordinates (30° elevation, looking down)
elev_rad = np.radians(self.elevation)
azim_rad = np.radians(self.azimuth)
cam_x = center[0] + camera_distance * np.cos(elev_rad) * np.cos(azim_rad)
cam_y = center[1] + camera_distance * np.cos(elev_rad) * np.sin(azim_rad)
cam_z = center[2] + camera_distance * np.sin(elev_rad)
⋮----
# Create renderer based on type
⋮----
renderer = TachyonRenderer()
⋮----
renderer = OSPRayRenderer()
⋮----
renderer = OpenGLRenderer()
⋮----
renderer = AnariRenderer()
⋮----
# Render to file
⋮----
# Clean up temporary file
⋮----
# Remove from scene
⋮----
"""Visualize a single structure."""
_ = title
_ = show
# Expand if needed
⋮----
output_file = 'output.png'
⋮----
"""Visualize trajectory as GIF animation."""
⋮----
output_file = str(Path(output_file).with_suffix('.gif'))
⋮----
n_frames = len(trajectory)
⋮----
titles = [f"Frame {i+1}/{n_frames}" for i in range(n_frames)]
⋮----
# Check if we need to expand structures
⋮----
trajectory = [make_supercell(atoms, P) for atoms in trajectory]
⋮----
temp_dir = tempfile.mkdtemp()
frame_files = []
⋮----
frame_file = os.path.join(temp_dir, f"frame_{i:04d}.png")
⋮----
images = [Image.open(f) for f in frame_files]
⋮----
def auto_detect_file_type(filename: str) -> Tuple[bool, int]
⋮----
"""Auto-detect if file is trajectory or single structure."""
⋮----
structures = read(filename, index=':')
⋮----
def main()
⋮----
"""Main CLI."""
parser = argparse.ArgumentParser(
⋮----
args = parser.parse_args()
⋮----
# Check input
⋮----
# Check OVITO availability
⋮----
# Parse supercell dimensions
⋮----
supercell = tuple(map(int, args.supercell.split(',')))
⋮----
# Detect file type
⋮----
# Create visualizer
visualizer = CatalystSurfaceVisualizer(
⋮----
# Output filename
⋮----
input_path = Path(args.input)
⋮----
# Info
⋮----
trajectory = read(args.input, index=':')
⋮----
atoms = read(args.input)
````

## File: core/viz/energy_diagram_plotter.py
````python
#!/usr/bin/env python
⋮----
"""
Energy Diagram Plotter - High-quality free energy diagram visualization
"""
⋮----
plt = None
rcParams = None
⋮----
def _ensure_matplotlib()
⋮----
plt = _plt
rcParams = _rc_params
⋮----
# Try to import catplot for smooth interpolation
⋮----
HAS_CATPLOT = True
⋮----
HAS_CATPLOT = False
⋮----
# =============================================================================
# Data Classes
⋮----
@dataclass
class EnergyStep
⋮----
"""Single step in energy diagram."""
label: str                          # Step label (e.g., 'H₂ + O₂', 'TS1')
energy: float                       # Energy value (kcal/mol, kJ/mol, eV, etc.)
is_ts: bool = False                 # Whether this is a transition state
color: Optional[str] = None         # Custom color for this step
line_style: Optional[str] = None    # Custom line style
⋮----
@dataclass
class DiagramStyle
⋮----
"""Styling configuration for energy diagram."""
# Colors
color_scheme: str = 'material'      # 'material', 'nature', 'science', 'custom'
primary_color: str = '#1976D2'      # Default blue
ts_color: str = '#D32F2F'           # Red for transition states
background_color: str = '#FFFFFF'   # White background
⋮----
# Line properties
line_width: float = 2.5
ts_line_width: float = 2.0
platform_length: float = 1.0        # Length of horizontal platforms
ts_peak_width: float = 0.8          # Width of TS peak
⋮----
# Shadow effects
use_shadow: bool = False            # Disabled by default
shadow_color: str = '#CCCCCC'
shadow_offset: float = 0.05         # Offset as fraction of y-range
⋮----
# Interpolation
interp_method: str = 'spline'       # 'spline' or 'quadratic'
interp_points: int = 100            # Points per segment
⋮----
# Labels
label_fontsize: int = 12
energy_fontsize: int = 10
show_energies: bool = True
energy_format: str = '.2f'
state_label_rotation: float = 0.0
⋮----
# Axes
xlabel: str = ''                    # No xlabel by default
ylabel: str = 'Free energy'
ylabel_units: str = 'kcal/mol'
grid: bool = False
⋮----
# Figure
figsize: Tuple[float, float] = (10, 6)
dpi: int = 300
⋮----
# Color Schemes
⋮----
COLOR_SCHEMES = {
⋮----
'primary': '#1976D2',      # Blue 700
'ts': '#D32F2F',           # Red 700
'secondary': '#388E3C',     # Green 700
'accent': '#F57C00',       # Orange 700
⋮----
'primary': '#0C4B8E',      # Nature blue
'ts': '#C3423F',           # Nature red
'secondary': '#7CB342',     # Nature green
'accent': '#FFB300',       # Nature yellow
⋮----
'primary': '#003C71',      # Science blue
'ts': '#B31B1B',           # Science red
'secondary': '#00A79D',     # Science teal
'accent': '#FF6F00',       # Science orange
⋮----
'primary': '#34495E',      # Dark gray-blue
'ts': '#E74C3C',           # Soft red
'secondary': '#16A085',     # Soft teal
'accent': '#F39C12',       # Soft orange
'background': '#ECF0F1',    # Light gray
'text': '#2C3E50',         # Dark blue-gray
⋮----
# Utility Functions
⋮----
def format_chemical_formula(formula: str) -> str
⋮----
"""
    Convert chemical formula to LaTeX format with subscripts and superscripts.

    Examples:
        'H2O' -> 'H$_2$O'
        'CO2' -> 'CO$_2$'
        'H+' -> 'H$^+$'
        'O2-' -> 'O$_2^-$'
        '1/2 O2' -> '$\\frac{1}{2}$ O$_2$'

    Args:
        formula: Chemical formula string

    Returns:
        LaTeX formatted string
    """
# Handle fractions
formula = re.sub(r'(\d+)/(\d+)', r'$\\frac{\1}{\2}$ ', formula)
⋮----
# Handle charges (+ or -)
formula = re.sub(r'(\d*)([+-])', r'$^{\1\2}$', formula)
⋮----
# Handle subscripts (numbers after letters)
formula = re.sub(r'([A-Za-z])(\d+)', r'\1$_{\2}$', formula)
⋮----
# Clean up multiple $ signs
formula = re.sub(r'\$+', '$', formula)
⋮----
"""
    Parse various input formats into list of EnergyStep objects.

    Args:
        energies: Energy values or dict/EnergyStep list
        labels: Optional labels for each step
        is_ts_list: Optional list indicating which steps are TS

    Returns:
        List of EnergyStep objects
    """
# Case 1: Already EnergyStep objects
⋮----
# Case 2: List of dictionaries
⋮----
# Case 3: List of floats
⋮----
steps = []
⋮----
# Auto-detect labels
⋮----
label = labels[i]
⋮----
label = f'State {i}'
⋮----
# Auto-detect TS
⋮----
is_ts = is_ts_list[i]
⋮----
# Heuristic: odd indices might be TS if energy is higher than neighbors
is_ts = False
⋮----
is_ts = True
⋮----
"""
    Generate smooth energy barrier curve using catplot's interpolation.

    Args:
        E_IS: Initial state energy
        E_TS: Transition state energy
        E_FS: Final state energy
        n_points: Number of interpolation points
        platform_length: Length of IS and FS platforms
        peak_width: Width of TS peak
        method: 'spline' or 'quadratic'

    Returns:
        (x, y) arrays for the barrier curve
    """
⋮----
# Use catplot's high-quality interpolation
⋮----
# Catplot fails for cases where TS is not highest (e.g., desorption)
# Fall back to scipy interpolation
⋮----
# Fallback: simple cubic spline (works for all energy profiles)
⋮----
# Define key points: IS platform -> barrier peak -> FS platform
x_key = np.array([
⋮----
0,                                          # Start of IS platform
platform_length,                            # End of IS platform / start of barrier
platform_length + peak_width / 2,          # TS peak position
platform_length + peak_width,              # End of barrier / start of FS platform
2 * platform_length + peak_width           # End of FS platform
⋮----
y_key = np.array([E_IS, E_IS, E_TS, E_FS, E_FS])
⋮----
# Interpolate
cs = CubicSpline(x_key, y_key, bc_type='clamped')
x = np.linspace(0, x_key[-1], n_points * 3)
y = cs(x)
⋮----
# Main Plotter Class
⋮----
class EnergyDiagramPlotter
⋮----
"""
    High-quality energy diagram plotter.

    Combines the best features of pMuTT and catplot to create publication-ready
    energy diagrams with minimal user input.

    Example:
        >>> plotter = EnergyDiagramPlotter()
        >>>
        >>> # Simple input: just energies and labels
        >>> energies = [0.0, 1.2, -0.5, 0.8, -2.5]
        >>> labels = ['H₂ + O₂', 'TS1', 'H₂O₂', 'TS2', '2 H₂O']
        >>>
        >>> plotter.plot(energies, labels)
        >>> plotter.save('energy_diagram.png', dpi=300)
    """
⋮----
def __init__(self, style: Optional[DiagramStyle] = None, apply_matplotlib_style: bool = True)
⋮----
"""
        Initialize plotter with optional custom style.

        Args:
            style: DiagramStyle object with custom styling
            apply_matplotlib_style: Whether to override global matplotlib rcParams
        """
⋮----
# Apply color scheme
⋮----
# Set matplotlib style
⋮----
def _apply_color_scheme(self)
⋮----
"""Apply predefined color scheme to style."""
⋮----
scheme = COLOR_SCHEMES[self.style.color_scheme]
⋮----
def _set_matplotlib_style(self)
⋮----
"""Set matplotlib rcParams for publication quality."""
⋮----
rc_params['text.usetex'] = False  # Can be set to True if LaTeX is available
⋮----
"""
        Create energy diagram.

        Args:
            energies: Energy values, dict list, or EnergyStep list
            labels: Optional labels for each step
            is_ts_list: Optional list indicating which steps are TS
            reference_index: Index of reference state (set to 0 energy)
            show: Whether to display the plot immediately
            ax: Optional matplotlib axis to draw into
        """
# Parse input
⋮----
# Set reference energy to zero
ref_energy = self.steps[reference_index].energy
⋮----
# Create figure
⋮----
# Plot energy profile
⋮----
# Format axes
⋮----
# Add labels
⋮----
def _plot_energy_profile(self)
⋮----
"""Plot the main energy profile with barriers using catplot's approach."""
x_offset = 0.0
x_positions = [None] * len(self.steps)
ts_positions = []
⋮----
# Color for barrierless connections
barrierless_color = '#999999'
barrierless_style = '--'
⋮----
# First pass: identify all barriers and handle the first state if it has no barrier before it
i = 0
⋮----
# Handle first state if it doesn't start with a barrier sequence
⋮----
first_step = self.steps[0]
# Check if first state is followed by a barrier
⋮----
# First state with no barrier after it - draw it as simple platform
⋮----
color = first_step.color if first_step.color else self.style.primary_color
⋮----
x_platform = np.array([x_offset, x_offset + self.style.platform_length])
y_platform = np.array([first_step.energy, first_step.energy])
⋮----
i = 1
⋮----
# Process all barriers
⋮----
step = self.steps[i]
⋮----
# Skip TS steps - they're handled as part of barriers
⋮----
# Check if next step is a TS (forming a barrier)
has_barrier = (i + 1 < len(self.steps) and
⋮----
# Draw barrier: IS (i) -> TS (i+1) -> FS (i+2)
ts_step = self.steps[i + 1]
fs_step = self.steps[i + 2]
⋮----
# Generate full reaction curve (IS platform + barrier + FS platform)
⋮----
# Check if this IS was already drawn as the FS of previous barrier
is_already_drawn = (x_positions[i] is not None)
⋮----
# Skip the IS platform part, only draw from barrier start
platform_end_idx = np.argmax(x_full >= self.style.platform_length)
x_plot = x_full[platform_end_idx:] - self.style.platform_length + x_offset
y_plot = y_full[platform_end_idx:]
⋮----
# Draw the full curve including IS platform
x_plot = x_full + x_offset
y_plot = y_full
# Store IS position
⋮----
# Find TS peak position
peak_idx = np.argmax(y_plot)
ts_x = x_plot[peak_idx]
ts_y = y_plot[peak_idx]
⋮----
# Store FS position (center of FS platform)
⋮----
# Determine color
color = step.color if step.color else self.style.primary_color
line_width = self.style.line_width
⋮----
# Plot the curve
⋮----
# Update x_offset to end of this curve
x_offset = x_plot[-1]
⋮----
# Move to FS (i+2) for next iteration
⋮----
# State without barrier after it
# Check if there's a next state to connect to
⋮----
# Draw current platform if not drawn yet
⋮----
y_platform = np.array([step.energy, step.energy])
⋮----
# Draw barrierless connection with dashed line
next_step = self.steps[i + 1]
gap = 0.3
x_next_start = x_offset + gap
⋮----
x_offset = x_next_start
⋮----
# Last state or followed by something else - just draw platform
⋮----
def _format_axes(self)
⋮----
"""Format axis labels, ticks, and grid."""
# X-axis
if self.style.xlabel:  # Only set if not empty
⋮----
self.ax.set_xticks([])  # Remove x-ticks for cleaner look
⋮----
# Y-axis
ylabel = self.style.ylabel
⋮----
# Grid
⋮----
# Spine styling
⋮----
def _add_labels(self)
⋮----
"""Add state labels and energy values."""
⋮----
# Add labels for non-TS states
⋮----
# Skip TS in label positions (they're part of barriers)
⋮----
# Format label (chemical formula with subscripts)
formatted_label = format_chemical_formula(step.label)
⋮----
# Add label below x-axis with appropriate spacing
y_range = self.ax.get_ylim()[1] - self.ax.get_ylim()[0]
y_label_pos = self.ax.get_ylim()[0] - y_range * 0.06
⋮----
# Add energy value above platform (closer to curve)
⋮----
energy_text = f'{step.energy:{self.style.energy_format}}'
⋮----
# Add labels for TS states (at peak positions)
⋮----
# Find the index of this TS in self.steps
ts_idx = None
⋮----
ts_idx = idx
⋮----
# Calculate activation barrier (energy difference from previous intermediate)
⋮----
# Find the previous non-TS state
prev_energy = None
⋮----
prev_energy = self.steps[j].energy
⋮----
# Display activation barrier instead of absolute energy
barrier = ts_step.energy - prev_energy
energy_text = f'{barrier:{self.style.energy_format}}'
⋮----
# Fallback to absolute energy if can't find previous state
energy_text = f'{ts_step.energy:{self.style.energy_format}}'
⋮----
# Fallback to absolute energy
⋮----
def save(self, filename: str, dpi: Optional[int] = None, **kwargs)
⋮----
"""
        Save figure to file.

        Args:
            filename: Output filename
            dpi: Resolution (overrides style.dpi)
            **kwargs: Additional arguments passed to plt.savefig
        """
⋮----
save_dpi = dpi or self.style.dpi
⋮----
def export_data(self, filename: str)
⋮----
"""
        Export energy data to CSV file.

        Args:
            filename: Output CSV filename
        """
⋮----
writer = csv.writer(f)
⋮----
# Convenience Functions
⋮----
"""
    Quick one-line energy diagram plotting.

    Args:
        energies: List of energy values
        labels: List of state labels
        is_ts_list: Optional list indicating TS states
        output: Output filename
        color_scheme: Color scheme name
        **style_kwargs: Additional style parameters

    Example:
        >>> quick_plot(
        ...     energies=[0, 1.2, -0.5, 0.8, -2.5],
        ...     labels=['H2 + O2', 'TS1', 'H2O2', 'TS2', '2 H2O'],
        ...     output='co_oxidation.png'
        ... )
    """
style = DiagramStyle(color_scheme=color_scheme, **style_kwargs)
plotter = EnergyDiagramPlotter(style=style)
⋮----
# Main entry point for testing
⋮----
# ==========================================================================
# Test 1: Simple CO oxidation pathway
⋮----
energies_1 = [0.0, 0.95, 0.15, 0.85, -2.35]
labels_1 = ['CO* + O*', 'TS1', 'CO-O*', 'TS2', 'CO₂ + 2*']
is_ts_1 = [False, True, False, True, False]
⋮----
plotter1 = EnergyDiagramPlotter(style=DiagramStyle(color_scheme='material'))
⋮----
# Test 2: H2O formation (using dict input)
⋮----
steps_2 = [
⋮----
style2 = DiagramStyle(
⋮----
plotter2 = EnergyDiagramPlotter(style=style2)
⋮----
# Test 3: Multiple color schemes comparison
⋮----
energies_3 = [0.0, 0.8, -0.5, 1.2, -1.8]
labels_3 = ['A', 'TS1', 'B', 'TS2', 'C']
⋮----
style = DiagramStyle(color_scheme=scheme, figsize=(8, 5))
⋮----
# Test 4: Quick plot function
⋮----
# Summary
````

## File: core/viz/visualization_manager.py
````python
"""
Unified Visualization Manager for Gas-Solid Catalysis Digital Twin

自动在以下关键步骤生成可视化：
1. SurFF 表面生成后 → 结构图片
2. AdsorbDiff 吸附位点后 → 结构图片
3. MC 表面重构模拟 → 轨迹 GIF
4. 反应路径中间体 → 完整反应 GIF
5. 能垒计算后 → 自由能图

Author: Claude
Date: 2026-01-21
"""
⋮----
# Import viz modules
⋮----
logger = logging.getLogger(__name__)
⋮----
class CatalysisVisualizationManager
⋮----
"""
    统一的可视化管理器

    在催化工作流的关键步骤自动生成高质量可视化：
    - 表面结构图
    - 吸附构型图
    - MC轨迹动画
    - 反应路径动画
    - 自由能图
    """
⋮----
"""
        Initialize visualization manager

        Parameters
        ----------
        output_dir : str
            输出目录
        quality : str
            渲染质量 ('low', 'medium', 'high', 'ultra')
        renderer : str
            渲染器 ('tachyon', 'ospray', 'opengl', 'anari')
        elevation : float
            视角仰角
        azimuth : float
            视角方位角
        color_scheme : str
            能量图配色方案 ('material', 'nature', 'science', 'elegant')
        enable_all : bool
            是否启用所有可视化（可以单独禁用某些步骤）
        """
⋮----
# 创建输出目录
⋮----
# 可视化开关
⋮----
# 初始化可视化器
⋮----
self.energy_plotter = None  # 延迟初始化
⋮----
# ========================================================================
# 1. SurFF 表面生成后可视化
⋮----
"""
        可视化 SurFF 生成的所有表面

        Parameters
        ----------
        surface_files : Dict[str, str]
            表面文件字典 {miller_index: filepath}
            例如: {'111': 'slab_111.vasp', '110': 'slab_110.vasp'}
        surface_energies : Dict[str, float], optional
            表面能 {miller_index: energy_eV_per_A2}

        Returns
        -------
        Dict[str, str]
            生成的图片路径 {miller_index: image_path}
        """
⋮----
surface_dir = self.viz_dir / "01_surfaces"
⋮----
output_paths = {}
⋮----
# 读取结构
slab = read(filepath)
⋮----
# 生成标题
title = f"Surface ({miller_idx})"
⋮----
energy = surface_energies[miller_idx]
⋮----
# 输出路径
output_path = surface_dir / f"surface_{miller_idx}.png"
⋮----
# 可视化
⋮----
# 2. AdsorbDiff 吸附位点后可视化
⋮----
"""
        可视化吸附构型

        Parameters
        ----------
        adsorbate_structures : Dict[str, str]
            吸附构型文件 {adsorbate_name: filepath}
            例如: {'CO': 'CO_adsorbed.vasp', 'O': 'O_adsorbed.vasp'}
        adsorption_energies : Dict[str, float], optional
            吸附能 {adsorbate_name: energy_eV}
        surface_name : str
            表面名称（用于文件命名）

        Returns
        -------
        Dict[str, str]
            生成的图片路径 {adsorbate_name: image_path}
        """
⋮----
ads_dir = self.viz_dir / "02_adsorption"
⋮----
structure = read(filepath)
⋮----
title = f"{ads_name} adsorbed"
⋮----
energy = adsorption_energies[ads_name]
⋮----
output_path = ads_dir / f"{surface_name}_{ads_name}_adsorbed.png"
⋮----
# 3. MC 表面重构轨迹可视化
⋮----
"""
        可视化 MC 表面重构轨迹为 GIF

        Parameters
        ----------
        trajectory_file : Union[str, List[Atoms]]
            轨迹文件路径(.traj)或直接传入的 Atoms 帧列表
        surface_name : str
            表面名称
        fps : int
            GIF 帧率
        max_frames : int, optional
            最大帧数（用于长轨迹采样）

        Returns
        -------
        str
            GIF 文件路径
        """
⋮----
mc_dir = self.viz_dir / "03_mc_reconstruction"
⋮----
frames = [atoms.copy() for atoms in trajectory_file if isinstance(atoms, Atoms)]
n_frames = len(frames)
⋮----
traj = Trajectory(trajectory_file)
n_frames = len(traj)
frames = [atoms for atoms in traj]
⋮----
indices = np.linspace(0, n_frames - 1, max_frames, dtype=int)
frames = [frames[i] for i in indices]
⋮----
titles = [f"MC Step {i+1}/{len(frames)}" for i in range(len(frames))]
output_path = mc_dir / f"{surface_name}_mc_reconstruction.gif"
⋮----
# 4. 反应路径完整过程可视化
⋮----
"""
        可视化完整反应路径为 GIF

        Parameters
        ----------
        pathway_structures : Dict[str, str]
            路径结构 {intermediate_name: filepath}
            例如: {'*CO': 'CO_ads.vasp', '*O': 'O_ads.vasp', '*CO2': 'CO2_ads.vasp'}
        sequence : List[str]
            反应顺序
            例如: ['*CO', '*O', '*CO2', '*']
        surface_name : str
            表面名称
        fps : int
            GIF 帧率
        hold_frames : int
            每个中间体停留帧数

        Returns
        -------
        str
            GIF 文件路径
        """
⋮----
pathway_dir = self.viz_dir / "04_reaction_pathway"
⋮----
# 读取所有结构
structures = {}
⋮----
# 构建帧序列（每个中间体重复几帧）
frames = []
titles = []
⋮----
output_path = pathway_dir / f"{surface_name}_reaction_pathway.gif"
⋮----
# 生成 GIF
⋮----
# 5. 自由能图可视化
⋮----
"""
        生成自由能图

        Parameters
        ----------
        energies : List[float]
            能量值列表 (eV)
        labels : List[str]
            状态标签
        is_ts_list : List[bool], optional
            过渡态标记
        barriers : Dict[str, float], optional
            能垒信息 {step_name: barrier_eV}
        filename : str
            输出文件名
        title : str, optional
            图表标题

        Returns
        -------
        str
            图片文件路径
        """
⋮----
energy_dir = self.viz_dir / "05_energy_diagram"
⋮----
# 创建样式
style = DiagramStyle(
⋮----
# 创建绘图器
plotter = EnergyDiagramPlotter(style=style)
⋮----
# 绘制
⋮----
# 添加标题（如果提供）
⋮----
# 保存
output_path = energy_dir / filename
⋮----
# 同时导出数据
csv_path = energy_dir / filename.replace('.png', '.csv')
⋮----
# 打印关键信息
⋮----
# 便捷方法：完整工作流可视化
⋮----
"""
        对完整的工作流结果进行可视化

        Parameters
        ----------
        workflow_result : WorkflowResult
            工作流的完整结果对象
        surface_name : str
            表面名称

        Returns
        -------
        Dict[str, str]
            所有生成的可视化文件路径
            {
                'surfaces': [path1, path2, ...],
                'adsorption': [path1, path2, ...],
                'mc_trajectory': path,
                'reaction_pathway': path,
                'energy_diagram': path
            }
        """
⋮----
viz_outputs = {
⋮----
# TODO: 根据实际的 workflow_result 结构调用各个可视化方法
# 这里提供模板
⋮----
# 1. 表面可视化
# if hasattr(workflow_result, 'surface_files'):
#     viz_outputs['surfaces'] = self.visualize_generated_surfaces(
#         workflow_result.surface_files,
#         workflow_result.surface_energies
#     )
⋮----
# 2. 吸附可视化
# if hasattr(workflow_result, 'adsorbate_structures'):
#     viz_outputs['adsorption'] = self.visualize_adsorption_sites(
#         workflow_result.adsorbate_structures,
#         workflow_result.adsorption_energies,
#         surface_name
⋮----
# 3. MC轨迹
# if hasattr(workflow_result, 'mc_trajectory_file'):
#     viz_outputs['mc_trajectory'] = self.visualize_mc_reconstruction(
#         workflow_result.mc_trajectory_file,
⋮----
# 4. 反应路径
# if hasattr(workflow_result, 'pathway_structures'):
#     viz_outputs['reaction_pathway'] = self.visualize_reaction_pathway(
#         workflow_result.pathway_structures,
#         workflow_result.intermediate_sequence,
⋮----
# 5. 能量图
# if hasattr(workflow_result, 'pathway_energies'):
#     viz_outputs['energy_diagram'] = self.visualize_energy_diagram(
#         workflow_result.pathway_energies,
#         workflow_result.intermediate_labels,
#         workflow_result.is_ts_list
⋮----
"""
        生成HTML摘要页面，展示所有可视化结果

        Parameters
        ----------
        viz_outputs : Dict[str, str]
            可视化输出路径
        output_file : str
            HTML文件名

        Returns
        -------
        str
            HTML文件路径
        """
⋮----
html_path = self.viz_dir / output_file
⋮----
html_content = f"""
⋮----
# 6. NEB轨迹可视化
⋮----
"""
        可视化NEB轨迹（插值路径或优化后的路径）

        Parameters
        ----------
        trajectory : List[Atoms]
            NEB轨迹帧列表
        reaction_name : str
            反应名称（如 "*CO_to_*CHO"）
        output_filename : str, optional
            输出文件名（默认自动生成）
        fps : int
            GIF帧率
        is_optimized : bool
            是否为优化后的路径（True）还是插值路径（False）

        Returns
        -------
        str
            GIF文件路径
        """
⋮----
neb_dir = self.viz_dir / "06_neb_trajectories"
⋮----
# 自动生成文件名
⋮----
prefix = "optimized" if is_optimized else "interpolated"
output_filename = f"{reaction_name}_{prefix}_path.gif"
⋮----
output_path = neb_dir / output_filename
⋮----
n_frames = len(trajectory)
⋮----
# 优化后的路径：提取最后n帧（收敛路径）
frame_titles = [
⋮----
# 插值路径
⋮----
# 生成GIF
⋮----
show_progress=False  # 不显示进度条，避免干扰日志
⋮----
"""
        从.traj文件读取并可视化NEB轨迹（提取收敛路径）

        Parameters
        ----------
        traj_file : str
            轨迹文件路径
        reaction_name : str, optional
            反应名称（从文件名推断）
        n_images : int
            提取最后n帧作为收敛路径
        output_filename : str, optional
            输出文件名
        fps : int
            GIF帧率

        Returns
        -------
        str
            GIF文件路径
        """
⋮----
# 读取轨迹
full_trajectory = ase_read(traj_file, index=':')
total_frames = len(full_trajectory)
⋮----
# 提取收敛路径（最后n帧）
⋮----
converged_path = full_trajectory[-n_images:]
⋮----
converged_path = full_trajectory
⋮----
# 推断反应名称
⋮----
reaction_name = Path(traj_file).stem.replace('_neb', '')
⋮----
# =============================================================================
# 测试代码
⋮----
def test_visualization_manager()
⋮----
"""测试可视化管理器"""
⋮----
# 创建测试目录
test_dir = "/tmp/viz_test"
⋮----
# 创建管理器
viz_manager = CatalysisVisualizationManager(
⋮----
# 1. 测试表面可视化
⋮----
slab_111 = fcc111('Pt', size=(4, 4, 4), vacuum=10.0)
slab_110 = fcc111('Pt', size=(4, 4, 4), vacuum=10.0)
⋮----
surface_files = {
surface_energies = {'111': 0.0928, '110': 0.1055}
⋮----
# 2. 测试吸附可视化
⋮----
slab_with_co = slab_111.copy()
co = molecule('CO')
⋮----
ads_structures = {'*CO': f"{test_dir}/CO_ads.vasp"}
ads_energies = {'*CO': -1.5}
⋮----
# 3. 测试能量图
⋮----
energies = [0.0, 1.2, 0.3, 0.9, -1.5]
labels = ['*CO + *O', 'TS1', '*CO-O', 'TS2', '*CO2 + *']
is_ts = [False, True, False, True, False]
````

"""Shared runtime utilities for CatDT CAMEL tools."""

import logging
import os
import pickle
from pathlib import Path
from typing import Any, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from camel_agents.gas_solid_digital_twin import GasSolidDigitalTwin
    from camel_agents.visualized_digital_twin import VisualizedGasSolidDigitalTwin
    from core.viz.visualization_manager import CatalysisVisualizationManager

CATDT_CORE_PATH = Path(__file__).resolve().parents[2] / "core"
DEPS_BASE_PATH = CATDT_CORE_PATH.parent / "deps"

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class CatDTToolRuntimeBase:
    def __init__(
        self,
        output_base_dir: str = "output/catdt_workflow",
        mc_energy_model: str = "CHGNet",
        potential_she: Optional[float] = None,
        ph: Optional[float] = None,
        mp_api_key: Optional[str] = None,
        adsorption_backend: str = "adsorbml",
        adsorbml_num_sites: int = 20,
        adsorbml_placement_mode: str = "random_site_heuristic_placement",
        adsorbml_interstitial_gap: float = 0.1,
        adsorbml_relax_steps: int = 200,
        adsorbml_relax_fmax: float = 0.05,
    ):
        self.output_base_dir = Path(output_base_dir)
        self.output_base_dir.mkdir(parents=True, exist_ok=True)
        self.dt_instance = None
        self.viz_manager_instance = None
        self._shared_fairchem_predictors: dict[str, Any] = {}

        # VSSR-MC energy model and electrochemical parameters
        self.mc_energy_model = mc_energy_model
        self.potential_she = potential_she
        self.ph = ph
        self.mp_api_key = mp_api_key

        # Agent2 adsorption backend selection
        self.adsorption_backend = adsorption_backend
        self.adsorbml_num_sites = adsorbml_num_sites
        self.adsorbml_placement_mode = adsorbml_placement_mode
        self.adsorbml_interstitial_gap = adsorbml_interstitial_gap
        self.adsorbml_relax_steps = adsorbml_relax_steps
        self.adsorbml_relax_fmax = adsorbml_relax_fmax

        # Set MP API key in environment if provided
        if mp_api_key:
            os.environ["MP_API_KEY"] = mp_api_key

    def _get_dt_instance(self, use_visualization: bool = False) -> Any:
        from camel_agents.gas_solid_digital_twin import GasSolidDigitalTwin
        from camel_agents.visualized_digital_twin import VisualizedGasSolidDigitalTwin

        if self.dt_instance is None or (use_visualization and not isinstance(self.dt_instance, VisualizedGasSolidDigitalTwin)):
            surff_root = str(DEPS_BASE_PATH / "SurFF")
            adsorbdiff_root = str(DEPS_BASE_PATH / "AdsorbDiff")
            surface_sampling_root = str(DEPS_BASE_PATH / "surface-sampling")
            fairchem_root = str(DEPS_BASE_PATH / "fairchem")
            fairchem_model_path = str(DEPS_BASE_PATH / "fairchem_models/uma-s-1p1.pt")

            # Map config mc_energy_model to VSSRMCPredictor model_type
            vssr_model = "CHGNetNFF"
            if self.mc_energy_model == "UMA":
                vssr_model = "UMA"

            shared_kwargs = dict(
                surff_root=surff_root,
                adsorbdiff_root=adsorbdiff_root,
                surface_sampling_root=surface_sampling_root,
                fairchem_root=fairchem_root,
                fairchem_model_path=fairchem_model_path,
                vssr_mc_model=vssr_model,
                potential_she=self.potential_she,
                ph=self.ph,
                use_gpu=True,
                adsorption_backend=self.adsorption_backend,
                adsorbml_num_sites=self.adsorbml_num_sites,
                adsorbml_placement_mode=self.adsorbml_placement_mode,
                adsorbml_interstitial_gap=self.adsorbml_interstitial_gap,
                adsorbml_relax_steps=self.adsorbml_relax_steps,
                adsorbml_relax_fmax=self.adsorbml_relax_fmax,
            )

            if use_visualization:
                self.dt_instance = VisualizedGasSolidDigitalTwin(
                    **shared_kwargs,
                    enable_llm_review=False,
                )
                logger.info("Initialized VisualizedGasSolidDigitalTwin (model=%s).", vssr_model)
            else:
                self.dt_instance = GasSolidDigitalTwin(
                    **shared_kwargs,
                    use_llm_controller=False,
                )
                logger.info("Initialized GasSolidDigitalTwin (model=%s).", vssr_model)
        return self.dt_instance

    def _get_viz_manager_instance(self, current_output_dir: Path) -> Any:
        from core.viz.visualization_manager import CatalysisVisualizationManager

        if self.viz_manager_instance is None or self.viz_manager_instance.output_dir != str(current_output_dir):
            self.viz_manager_instance = CatalysisVisualizationManager(
                output_dir=str(current_output_dir),
                quality='high',
                renderer='tachyon'
            )
            logger.info(f"Initialized CatalysisVisualizationManager for {current_output_dir}.")
        return self.viz_manager_instance

    def _make_pickle_safe(self, obj: Any, depth: int = 0) -> Any:
        if depth > 5:
            return repr(obj)

        if isinstance(obj, (str, int, float, bool, type(None))):
            return obj

        if isinstance(obj, dict):
            return {
                str(k): self._make_pickle_safe(v, depth + 1)
                for k, v in list(obj.items())[:500]
            }

        if isinstance(obj, (list, tuple, set)):
            return [self._make_pickle_safe(v, depth + 1) for v in list(obj)[:500]]

        if hasattr(obj, "tolist"):
            try:
                return self._make_pickle_safe(obj.tolist(), depth + 1)
            except Exception:
                pass

        if hasattr(obj, "item"):
            try:
                return obj.item()
            except Exception:
                pass

        if hasattr(obj, "__dict__"):
            safe = {"__type__": type(obj).__name__}
            for key, value in vars(obj).items():
                if str(key).startswith("_"):
                    continue
                safe[str(key)] = self._make_pickle_safe(value, depth + 1)
            return safe

        return repr(obj)

    def _save_result_to_pickle(self, result: Any, path: Path) -> str:
        """Helper to save complex objects to pickle file."""
        with open(path, 'wb') as f:
            try:
                pickle.dump(result, f)
            except Exception as exc:
                logger.warning("Pickle dump failed for %s, saving sanitized fallback: %s", path, exc)
                fallback = {
                    "_pickle_fallback": True,
                    "error": str(exc),
                    "result_type": type(result).__name__,
                    "result": self._make_pickle_safe(result),
                }
                pickle.dump(fallback, f)
        return str(path)

    def _load_result_from_pickle(self, path: str) -> Any:
        """Helper to load complex objects from pickle file."""
        with open(path, 'rb') as f:
            return pickle.load(f)

    def get_shared_fairchem_predictor(
        self,
        cache_key: str,
        model_name: str = "uma-s-1p1",
        use_gpu: bool = True,
        device: str = "cuda",
        work_subdir: str = "_shared_fairchem",
        keep_files: bool = False,
        verbose: bool = False,
    ) -> Any:
        predictor = self._shared_fairchem_predictors.get(cache_key)
        if predictor is not None:
            return predictor

        from core.pathway.fairchem_predictor import FairchemPredictor

        predictor_work_dir = self.output_base_dir / work_subdir / cache_key
        predictor_work_dir.mkdir(parents=True, exist_ok=True)
        predictor = FairchemPredictor(
            fairchem_root=str(DEPS_BASE_PATH / "fairchem"),
            model_name=model_name,
            use_gpu=use_gpu,
            device=device,
            work_dir=str(predictor_work_dir),
            keep_files=keep_files,
            verbose=verbose,
        )
        predictor._load_model()
        self._shared_fairchem_predictors[cache_key] = predictor
        logger.info("Initialized shared FairchemPredictor cache key=%s", cache_key)
        return predictor

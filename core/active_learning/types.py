from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np


@dataclass
class CandidateScore:
    structure_id: str
    uncertainty_force: float
    uncertainty_energy: float
    uncertainty_potential: float = 0.0
    descriptor: Optional[np.ndarray] = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def uncertainty_total(self) -> float:
        return float(
            self.uncertainty_force
            + self.uncertainty_energy
            + self.uncertainty_potential
        )

"""An experiment = (code, noise, decoder, cutoff, metrics), stored as plain data so results carry their config."""

import copy
import itertools
from dataclasses import asdict, dataclass, field, replace

from cvqec.circuits.teleportation import TeleportationParams


@dataclass
class Experiment:
    code: str = "gkp"
    code_params: dict = field(default_factory=lambda: {"lattice": "square", "delta": 0.3})
    decoder: str = "auto"
    decoder_params: dict = field(default_factory=dict)
    noise: TeleportationParams = field(default_factory=lambda: TeleportationParams.from_db(10.0))
    N: int = 140
    metrics: tuple[str, ...] = ("average_fidelity", "entanglement_fidelity", "pauli_errors", "mean_photon")
    label: str = ""

    def flat(self) -> dict:
        """One flat dict (for a CSV row): code.lattice, noise.eta, ..."""
        out = {"label": self.label or self.default_label(), "code": self.code, "decoder": self.decoder, "N": self.N}
        out.update({f"code.{k}": v for k, v in self.code_params.items()})
        out.update({f"decoder.{k}": v for k, v in self.decoder_params.items()})
        out.update({f"noise.{k}": v for k, v in asdict(self.noise).items()})
        return out

    def default_label(self) -> str:
        extra = self.code_params.get("lattice", "")
        return f"{self.code}{'_' + extra if extra else ''}/{self.decoder}"

    def with_(self, path: str, value) -> "Experiment":
        """Copy with one field changed; `path` is 'N', 'decoder', 'noise.eta', 'code_params.delta', ..."""
        exp = copy.deepcopy(self)
        if "." not in path:
            return replace(exp, **{path: value})
        head, key = path.split(".", 1)
        if head == "noise":
            exp.noise = replace(exp.noise, **{key: value})
        else:
            getattr(exp, head)[key] = value
        return exp


def grid(base: Experiment, **axes) -> list[Experiment]:
    """Cartesian product over axes, e.g. grid(base, **{"noise.eta": [0.9, 1.0], "decoder": ["petz", "auto"]})."""
    names = list(axes)
    exps = []
    for values in itertools.product(*(axes[n] for n in names)):
        exp = base
        for n, v in zip(names, values):
            exp = exp.with_(n, v)
        exps.append(exp)
    return exps

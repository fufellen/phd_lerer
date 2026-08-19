"""Оптические постоянные для трёхмерных сцен, с указанием источника.

Каждая константа сопровождается ссылкой на работу или на файл в этом
репозитории, откуда она взята. Там, где для одного материала в разных задачах
приняты разные наборы (золото у Park и золото в статье про PCM), сохранены оба,
и сцена берёт тот, с которым сверяется её результат. Смешивать наборы нельзя:
расхождение по Re(n_eff) плёнки золота между ними составляет около 1e-3.

Соглашение: eps = (n + i k)^2 при k > 0, поэтому у поглощающей среды Im(eps) > 0
и у затухающей моды Im(n_eff) > 0. То же соглашение принято в lrspp_coupling и
в edp/metal_strip_w800.
"""

from __future__ import annotations

from dataclasses import dataclass


def eps_from_nk(n: float, k: float = 0.0) -> complex:
    """eps = (n + i k)^2."""
    return complex(n, k) ** 2


def k0_from_lambda(lambda_um: float) -> float:
    """Волновое число свободного пространства, мкм^-1."""
    from math import pi

    return 2.0 * pi / lambda_um


# --------------------------------------------------------------- металлы

# Набор Park: Palik, приведён в Park et al. 2011, рис. 1. Им пользуется
# lrspp_coupling, поэтому сцены 1-3 берут его.
AU_PARK_1550 = eps_from_nk(0.550, 11.4912)
AU_PARK_1310 = eps_from_nk(0.42, 9.1)

# Набор статьи про PCM LR-DLSPP: воспроизводит референсные моды Si/Au и Si/Au/Si
# из data/au_thickness_sweep_reference.csv, поэтому сцена 4 берёт его.
AU_ARTICLE_1550 = eps_from_nk(0.6389, 11.1748)

# Серебро для видимого диапазона: Johnson & Christy 1972.
AG_633 = eps_from_nk(0.135, 3.9880)
AU_633 = eps_from_nk(0.1834, 3.4332)


# ------------------------------------------------------------ диэлектрики

N_ZPU450 = 1.450                     # Park 2011, Lee 2019: обкладка и центр
ZPU450 = eps_from_nk(N_ZPU450)
N_ZPU13_430 = 1.43                   # Park 2009: обкладка
N_ZPU13_440 = 1.44                   # Park 2009: сердцевина
ZPU13_430 = eps_from_nk(N_ZPU13_430)
ZPU13_440 = eps_from_nk(N_ZPU13_440)

SI_1550 = eps_from_nk(3.478)         # c-Si при 1550 нм
SIO2_1550 = eps_from_nk(1.444)       # SiO2 при 1550 нм
AIR = complex(1.0, 0.0)


# ------------------------------------------------------ материалы с фазовым переходом


@dataclass(frozen=True)
class PhaseChangeMaterial:
    """Материал с фазовым переходом в двух состояниях при 1550 нм."""

    name: str
    n_amorphous: float
    k_amorphous: float
    n_crystalline: float
    k_crystalline: float
    source: str

    @property
    def eps_amorphous(self) -> complex:
        return eps_from_nk(self.n_amorphous, self.k_amorphous)

    @property
    def eps_crystalline(self) -> complex:
        return eps_from_nk(self.n_crystalline, self.k_crystalline)

    @property
    def delta_n_material(self) -> float:
        """Изменение показателя самого материала, не моды."""
        return self.n_crystalline - self.n_amorphous

    def eps(self, state: str) -> complex:
        if state == "amorphous":
            return self.eps_amorphous
        if state == "crystalline":
            return self.eps_crystalline
        raise ValueError(f"состояние должно быть amorphous или crystalline, получено {state!r}")


# Значения и источники повторяют data/pcm_optical_constants.csv статьи про
# ЭДП PCM LR-DLSPP, чтобы сцена сверялась с уже проверенными там числами.
PCM = {
    "GSST": PhaseChangeMaterial(
        "Ge2Sb2Se4Te1 (GSST)",
        3.47, 0.0002,
        5.5, 0.42,
        "Gosciniak 2022 doi:10.1063/5.0082094, таблица [17]",
    ),
    "Sb2S3": PhaseChangeMaterial(
        "Sb2S3",
        2.712, 1e-5,
        3.308, 1e-5,
        "Delaney 2020 doi:10.1002/adfm.202002447; k < 1e-5 указано как оценка сверху",
    ),
    "Sb2Se3": PhaseChangeMaterial(
        "Sb2Se3",
        3.285, 1e-5,
        4.055, 1e-5,
        "Delaney 2020 doi:10.1002/adfm.202002447; k < 1e-5 указано как оценка сверху",
    ),
}

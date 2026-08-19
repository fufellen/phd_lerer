"""Аналитические расчёты плазмонных слоистых волноводов.

Ядро - комплексный решатель TM-мод произвольной слоистой структуры методом
матрицы переноса (tmm), поверх него метод эффективного показателя для полоски
конечной ширины (eim) и модель связанных волноводов (cmt).

Используется для воспроизведения результатов статей:
- Park et al., Opt. Commun. 282, 4513 (2009) - вертикальный направленный
  ответвитель между LR-SPP и диэлектрическим волноводом;
- Park et al., Opt. Express 19, 21605 (2011) - плазмонный преобразователь
  размера моды ДМД -> ДМДМД.
"""

from .tmm import (
    Stack,
    find_modes,
    mode_fields,
    mode_width_1e,
    orthogonality_residual,
    overlap_power,
    power_flux,
    propagation_length_um,
    propagation_loss_db_per_cm,
    solve_mode,
    tm_residual,
    tm_residual_3layer,
)
from .eim import solve_strip, strip_cutoff_width
from .cmt import CoupledPair, CouplerModel, extract_kappa, optimal_length_um
from . import materials

__all__ = [
    "Stack",
    "find_modes",
    "mode_fields",
    "mode_width_1e",
    "orthogonality_residual",
    "overlap_power",
    "power_flux",
    "propagation_length_um",
    "propagation_loss_db_per_cm",
    "solve_mode",
    "tm_residual",
    "tm_residual_3layer",
    "solve_strip",
    "strip_cutoff_width",
    "CoupledPair",
    "CouplerModel",
    "extract_kappa",
    "optimal_length_um",
    "materials",
]

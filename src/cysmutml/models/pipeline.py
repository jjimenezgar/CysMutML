"""Scikit-learn model pipelines."""

from __future__ import annotations

from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


def preprocessor(
    numeric: list[str], categorical: list[str], scale_numeric: bool
) -> ColumnTransformer:
    numeric_steps = [("imputer", SimpleImputer(strategy="median"))]
    if scale_numeric:
        numeric_steps.append(("scaler", StandardScaler()))
    categorical_steps = [
        ("imputer", SimpleImputer(strategy="most_frequent")),
        ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
    ]
    return ColumnTransformer(
        [
            ("numeric", Pipeline(numeric_steps), numeric),
            ("categorical", Pipeline(categorical_steps), categorical),
        ],
        remainder="drop",
    )


def make_regressors(
    numeric: list[str], categorical: list[str], random_seed: int
) -> dict[str, Pipeline]:
    return {
        "dummy_mean": Pipeline(
            [("preprocess", preprocessor(numeric, categorical, False)), ("model", DummyRegressor())]
        ),
        "ridge": Pipeline(
            [("preprocess", preprocessor(numeric, categorical, True)), ("model", Ridge(alpha=1.0))]
        ),
        "random_forest": Pipeline(
            [
                ("preprocess", preprocessor(numeric, categorical, False)),
                (
                    "model",
                    RandomForestRegressor(
                        n_estimators=100, random_state=random_seed, min_samples_leaf=2
                    ),
                ),
            ]
        ),
        "hist_gradient_boosting": Pipeline(
            [
                ("preprocess", preprocessor(numeric, categorical, False)),
                ("model", HistGradientBoostingRegressor(random_state=random_seed, max_iter=100)),
            ]
        ),
    }

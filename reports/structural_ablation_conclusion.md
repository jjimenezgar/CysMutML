# Structural Ablation Conclusion

Status: PARTIALLY EXECUTED.

This run used a balanced-small structural candidate subset, not the full FireProtDB structural candidate set. Therefore the conclusion is preliminary.

## Main Question

Does explicit structural context improve prediction over physicochemical descriptors?

Current evidence on the executed subset: **no detectable improvement**. Structural features worsened mean MAE for Ridge and HistGradientBoosting on these rows.

Ridge paired delta MAE, structure minus physchem:

model  fold  delta_mae_structure_minus_physchem  delta_rmse_structure_minus_physchem  delta_r2_structure_minus_physchem
ridge     1                            0.290988                             0.181866                          -0.338198
ridge     2                            0.878527                             1.189311                          -2.516592
ridge     3                           -0.080354                             0.051014                          -0.035383

HistGradientBoosting paired delta MAE, structure minus physchem:

                 model  fold  delta_mae_structure_minus_physchem  delta_rmse_structure_minus_physchem  delta_r2_structure_minus_physchem
hist_gradient_boosting     1                            0.556824                             0.416536                          -0.548671
hist_gradient_boosting     2                            0.015478                             0.296315                          -0.463352
hist_gradient_boosting     3                           -0.010500                             0.131918                          -0.088099

Interpretation: positive delta MAE means structural features were worse. The subset is small (114 rows, 60 proteins), so this should not be generalized as a final negative result. It is evidence that the current structural subset/mapping/features are not yet ready to replace the deployed physicochemical Ridge model.

Decision: do not promote a structural model artifact.

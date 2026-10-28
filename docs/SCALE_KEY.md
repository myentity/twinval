# Scale and sub-scale key

Generated from `src/twinval/registry/instruments.yaml` (`python tools/make_scale_key.py`).
Column *items* lists the item codes **in the order of the `it1 ... itK` columns** of the paired data files.
*Redundant* scales (facets nested in a domain, total scores) are analysed at the scale level but are excluded from the
profile, ESEM and item-level analyses, where they would count the same items twice.

### 148 - Big Five Inventory-2 (BFI-2) - Soto & John (2017)

Response scale 1-5; item prefix `BFI2`.

| scale_id | name | level | scoring | items (it1...itK) | reverse-keyed items | in profile / ESEM / item level |
|---|---|---|---|---|---|---|
| `148_Extraversion` | Extraversion | domain | mean | BFI2_1, BFI2_6, BFI2_11, BFI2_16, BFI2_21, BFI2_26, BFI2_31, BFI2_36, BFI2_41, BFI2_46, BFI2_51, BFI2_56 | BFI2_11, BFI2_16, BFI2_26, BFI2_31, BFI2_36, BFI2_51 | yes |
| `148_Sociability` | Sociability | facet | mean | BFI2_1, BFI2_16, BFI2_31, BFI2_46 | BFI2_16, BFI2_31 | no (redundant: nested in `148_Extraversion`) |
| `148_Persistence` | Assertiveness | facet | mean | BFI2_6, BFI2_21, BFI2_36, BFI2_51 | BFI2_36, BFI2_51 | no (redundant: nested in `148_Extraversion`) |
| `148_EnergyLevel` | Energy Level | facet | mean | BFI2_11, BFI2_26, BFI2_41, BFI2_56 | BFI2_11, BFI2_26 | no (redundant: nested in `148_Extraversion`) |
| `148_Agreeableness` | Agreeableness | domain | mean | BFI2_2, BFI2_7, BFI2_12, BFI2_17, BFI2_22, BFI2_27, BFI2_32, BFI2_37, BFI2_42, BFI2_47, BFI2_52, BFI2_57 | BFI2_12, BFI2_17, BFI2_22, BFI2_37, BFI2_42, BFI2_47 | yes |
| `148_Compassion` | Compassion | facet | mean | BFI2_2, BFI2_17, BFI2_32, BFI2_47 | BFI2_17, BFI2_47 | no (redundant: nested in `148_Agreeableness`) |
| `148_Respectfulness` | Respectfulness | facet | mean | BFI2_7, BFI2_22, BFI2_37, BFI2_52 | BFI2_22, BFI2_37 | no (redundant: nested in `148_Agreeableness`) |
| `148_Trust` | Trust | facet | mean | BFI2_12, BFI2_27, BFI2_42, BFI2_57 | BFI2_12, BFI2_42 | no (redundant: nested in `148_Agreeableness`) |
| `148_Conscientiousness` | Conscientiousness | domain | mean | BFI2_3, BFI2_8, BFI2_13, BFI2_18, BFI2_23, BFI2_28, BFI2_33, BFI2_38, BFI2_43, BFI2_48, BFI2_53, BFI2_58, BFI2_61 | BFI2_3, BFI2_8, BFI2_23, BFI2_28, BFI2_48, BFI2_58, BFI2_61 | yes |
| `148_Organization` | Organization | facet | mean | BFI2_3, BFI2_18, BFI2_33, BFI2_48 | BFI2_3, BFI2_48 | no (redundant: nested in `148_Conscientiousness`) |
| `148_Productiveness` | Productiveness | facet | mean | BFI2_8, BFI2_23, BFI2_38, BFI2_53 | BFI2_8, BFI2_23 | no (redundant: nested in `148_Conscientiousness`) |
| `148_Responsibility` | Responsibility | facet | mean | BFI2_13, BFI2_28, BFI2_43, BFI2_58, BFI2_61 | BFI2_28, BFI2_58, BFI2_61 | no (redundant: nested in `148_Conscientiousness`) |
| `148_NegativeEmotionality` | Neuroticism | domain | mean | BFI2_4, BFI2_9, BFI2_14, BFI2_19, BFI2_24, BFI2_29, BFI2_34, BFI2_39, BFI2_44, BFI2_49, BFI2_54, BFI2_59 | BFI2_4, BFI2_9, BFI2_24, BFI2_29, BFI2_44, BFI2_49 | yes |
| `148_Anxiety` | Anxiety | facet | mean | BFI2_4, BFI2_19, BFI2_34, BFI2_49 | BFI2_4, BFI2_49 | no (redundant: nested in `148_NegativeEmotionality`) |
| `148_Depression` | Depression | facet | mean | BFI2_9, BFI2_24, BFI2_39, BFI2_54 | BFI2_9, BFI2_24 | no (redundant: nested in `148_NegativeEmotionality`) |
| `148_EmotionalVolatility` | Emotional Volatility | facet | mean | BFI2_14, BFI2_29, BFI2_44, BFI2_59 | BFI2_29, BFI2_44 | no (redundant: nested in `148_NegativeEmotionality`) |
| `148_OpenMindedness` | Openness to Experience | domain | mean | BFI2_5, BFI2_10, BFI2_15, BFI2_20, BFI2_25, BFI2_30, BFI2_35, BFI2_40, BFI2_45, BFI2_50, BFI2_55, BFI2_60 | BFI2_5, BFI2_25, BFI2_30, BFI2_45, BFI2_50, BFI2_55 | yes |
| `148_AestheticSensitivity` | Aesthetic Sensitivity | facet | mean | BFI2_5, BFI2_20, BFI2_35, BFI2_50 | BFI2_5, BFI2_50 | no (redundant: nested in `148_OpenMindedness`) |
| `148_IntellectualCuriosity` | Intellectual Curiosity | facet | mean | BFI2_10, BFI2_25, BFI2_40, BFI2_55 | BFI2_25, BFI2_55 | no (redundant: nested in `148_OpenMindedness`) |
| `148_CreativeImagination` | Creative Imagination | facet | mean | BFI2_15, BFI2_30, BFI2_45, BFI2_60 | BFI2_30, BFI2_45 | no (redundant: nested in `148_OpenMindedness`) |

### 161 - Mental Health Continuum - Keyes (2009), Mental Health Continuum - Short Form (MHC-SF)

Response scale 1-6; item prefix `SPH`.

| scale_id | name | level | scoring | items (it1...itK) | reverse-keyed items | in profile / ESEM / item level |
|---|---|---|---|---|---|---|
| `161_EmotionalWellBeing` | Emotional Well-Being | subscale | sum | SPH_1, SPH_2, SPH_3 | - | yes |
| `161_SocialWellBeing` | Social Well-Being | subscale | sum | SPH_4, SPH_5, SPH_6, SPH_7, SPH_8 | - | yes |
| `161_PsychologicalWellBeing` | Psychological Well-Being | subscale | sum | SPH_9, SPH_10, SPH_11, SPH_12, SPH_13, SPH_14 | - | yes |
| `161_OverallScore` | Overall Score | total | sum | SPH_1, SPH_2, SPH_3, SPH_4, SPH_5, SPH_6, SPH_7, SPH_8, SPH_9, SPH_10, SPH_11, SPH_12, SPH_13, SPH_14 | - | no (redundant: total score) |

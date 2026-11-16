# RTMB spatially varying coefficients match frozen TMB values

    list(key = "effects svc spatial on", fn = 62.63126559, gr = c(-0.07050529229, 
    0.1187005359, -1.820294857, -1.544744677, -4.154613052, 50.36763701
    ), joint = list(fn = -51.43440223, gr = c(90, 57.78661034, 2038.074275, 
    -214.3679526)))

---

    list(key = "effects svc projection spatial on", report = list(
        proj_eta = c(0.4701361899, 0.5888482663, 1.102458728, 0.1190595283, 
        0.4184354424)))

---

    list(key = "effects svc spatial off", fn = 62.29290148, gr = c(-0.1040923008, 
    0.1213528712, -1.572107746, -1.546546013, -4.078192346, 50.58215319
    ), joint = list(fn = -84.88326059, gr = c(90, 7.267521997, 1951.35611, 
    -42.09067639)))

---

    list(key = "effects svc projection spatial off", report = list(
        proj_eta = c(0.4669428246, 0.5816179315, 1.100967635, 0.1204333573, 
        0.4147012404)))

# RTMB time-varying coefficients match frozen TMB values

    list(key = "effects time-varying rw", fn = 84.79285646, gr = c(-5.253906327, 
    6.751765808, 1.239420413, 1.519700543), joint = list(fn = 84.65034867, 
        gr = c(-4.384298813, 4.307941036, 1.964444444, 1.964444444, 
        -5.481049037, 0.1916244663, -0.05091025857, 3.515042551, 
        5.756703954, 7.252042817)))

---

    list(key = "effects time-varying projection rw", report = list(
        proj_eta = c(0.5133886231, -0.06534569719, 0.2346754779, 
        0.1665952721, -0.09602131447), proj_rw_i = c(0.5133886231, 
        -0.06534569719, 0.2346754779, 0.1665952721, -0.09602131447
        )))

---

    list(key = "effects time-varying rw0", fn = 85.88874492, gr = c(-5.330475736, 
    6.908031229, 2.084281649, 2.290898093), joint = list(fn = 85.52435226, 
        gr = c(-4.384298813, 4.307941036, 2.853333333, 2.96, -6.036604593, 
        0.1916244663, -0.05091025857, 3.626153663, 5.756703954, 7.252042817
        )))

---

    list(key = "effects time-varying projection rw0", report = list(
        proj_eta = c(0.4324904504, -0.05551477533, 0.1980703065, 
        0.1429118403, -0.07610749321), proj_rw_i = c(0.4324904504, 
        -0.05551477533, 0.1980703065, 0.1429118403, -0.07610749321
        )))

---

    list(key = "effects time-varying ar1", fn = 85.48491834, gr = c(-5.649137097, 
    6.759101541, 2.361051171, 2.371809613, -0.07404281257, -0.185587409
    ), joint = list(fn = 85.48216909, gr = c(12, 16.40477784, 38.29661982, 
    -17.30319893)))

---

    list(key = "effects time-varying projection ar1", report = list(
        proj_eta = c(0.4834968951, -0.05787364498, 0.2180666232, 
        0.1366867336, -0.1280686522), proj_rw_i = c(0.4834968951, 
        -0.05787364498, 0.2180666232, 0.1366867336, -0.1280686522
        )))

# RTMB correlated IID effects match frozen TMB values

    list(key = "effects correlated iid", fn = 55.36907918, gr = c(13, 
    23.61584034, 27.46874596, 2.246343069), joint = list(fn = 55.29087899, 
        gr = c(33, 30.99679043, 74.5646403, 2.44709555)))

---

    list(key = "effects correlated iid projection", report = list(
        proj_eta = c(0.0680176978, 0.2388946588, 0.04916816486, -0.4939915662, 
        -0.1904276351), proj_iid_re_i = c(0.0680176978, 0.2388946588, 
        0.04916816486, -0.4939915662, -0.1904276351)))

# RTMB thresholds and their priors match frozen TMB values

    list(key = "effects threshold breakpt", fn = 102.4927009, gr = c(-19.38125007, 
    -30.50702373, 39.11617089, 4.549635561))

---

    list(key = "effects threshold projection breakpt", report = list(
        proj_eta = c(-1.4, -1.352542373, -1.305084746, -1.257627119, 
        -1.210169492)))

---

    list(key = "effects threshold logistic", fn = 110.1109007, gr = c(35.74909834, 
    -42.58925183, -9.930633383, -2.794574241, 35.55380876))

---

    list(key = "effects threshold projection logistic", report = list(
        proj_eta = c(0.001100887236, 0.001412533898, 0.001812270262, 
        0.002324909556, 0.002982198886)))

# RTMB Matérn and transformed-parameter priors match frozen TMB values

    list(key = "effects matern priors bayesian FALSE", fn = 107.0958961, 
        gr = c(0.6336021651, -4.96358316, -1.31134609, -2.120670972, 
        -7.888744328, -64.49656288, -0.05325696016), joint = list(
            fn = -276.9851022, gr = c(187, -170.7554801, 460.4733128, 
            166.3334543)))

---

    list(key = "effects matern priors bayesian TRUE", fn = 110.4767579, 
        gr = c(0.6336021651, -4.96358316, -0.3113460902, -1.120670972, 
        -4.888744328, -65.49656288, 0.1441183601), joint = list(fn = -273.6042405, 
            gr = c(187, -166.5581048, 456.2759375, 163.2500874)))

# RTMB REML and profiled fixed effects fit

    list(key = "effects reml", fn = 50.1616664, gr = c(-2.193218834, 
    -1.964175287, 51.1398946), joint = list(fn = 165.590133, gr = c(44, 
    117.5102945, 3133.592464, -437.2075673)))

# RTMB dispersion formulas match frozen TMB values

    list(key = "effects dispformula gaussian", fn = 91.24946952, 
        gr = c(6.685202894, 2.449871801, -31.04387541, 34.38151037
        ))

---

    list(key = "effects dispformula nbinom1", fn = 104.6501884, gr = c(-16.4503799, 
    -11.82774506, -0.5070513578, 3.563190322))

---

    list(key = "effects dispformula nbinom2", fn = 108.6066873, gr = c(-17.76443947, 
    -8.307789927, -5.768907095, 1.43356231))

# RTMB gamma prior on time-varying SD matches frozen TMB values

    list(key = "effects gamma prior sigma_V", fn = 100.1524239, gr = c(-0.6082108695, 
    -2.002407024, -3.549980498, -14.13830757, -0.2606140492, 7.235466372, 
    1.605432446), joint = list(fn = -170.568854, gr = c(154, -172.8191096, 
    317.1972951, -78.33308213)))

# RTMB combined model matches every TMB report and sdreport row

    list(key = "effects combined", fn = 121.0338575, gr = c(13, 1.056174968, 
    38.18746573, 10.10094421), joint = list(fn = -92.52685617, gr = c(364, 
    252.1286065, 10323.66133, 169.8234124)))

---

    list(key = "effects combined projection", report = list(proj_eta = c(0.7749413968, 
    -0.6046358804, -0.5940711363)))


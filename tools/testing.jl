using Clapeyron

model = SRK(["limonene", "octanoic acid", "ethanol", "acetic acid"],
    userlocations = ["./database/"],
)

model.alpha.params.acentricfactor
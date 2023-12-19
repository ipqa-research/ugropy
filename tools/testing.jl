using Clapeyron

model = PSRK(["limonene"],
    userlocations=["tools/database/"],
    group_userlocations = ["tools"],
)

model.mixing.activity
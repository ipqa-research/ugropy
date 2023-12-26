using Clapeyron

model = PSRK(["limonene", "octanoic acid"],
    userlocations=["tools/database/"],
    group_userlocations = ["tools/database/PSRK/"],
    translation_userlocations = ["tools/database/"]
)

model.mixing.activity

model = ogUNIFAC(
    ["limonene", "octanoic acid"],
    group_userlocations = ["tools/database/ogUNIFAC/"],
    )
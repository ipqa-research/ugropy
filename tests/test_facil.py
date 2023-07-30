import ugropy as ug


def test_facil():
    subs = ug.Substance("hexane")

    print(subs.unifac_groups)
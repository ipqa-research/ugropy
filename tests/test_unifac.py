from cases import TCase

import pytest

from ugropy import unifac


@pytest.mark.dependency(name="unifac")
@pytest.mark.unifac
class TestUNIFAC(TCase):
    # Store al the groups detected in all the cases here:
    tested_groups = set()

    def asserts(self, case, solver):
        if case.unifac_result is None:
            pytest.skip(
                f"No UNIFAC result defined for {case.identifier}, "
                f"{case.cases_module}"
            )

        result = unifac.get_groups(
            identifier=case.identifier,
            identifier_type=case.identifier_type,
            solver=solver,
            search_multiple_solutions=True,
        )

        if len(result) > 1:
            for r in result:
                comp = r.subgroups in case.unifac_result

                if not comp:
                    message = (
                        "\n"
                        f"Case: {case.identifier}\n"
                        f"Test module: {case.cases_module}\n"
                        f"Expected: {case.unifac_result}\n"
                        f"Obtained: {[r.subgroups for r in result]}"
                    )

                    assert False, message

                self.tested_groups.update(r.subgroups.keys())
        else:
            comp = result[0].subgroups == case.unifac_result

            self.tested_groups.update(result[0].subgroups.keys())

            if not comp:
                message = (
                    "\n"
                    f"Case: {case.identifier}\n"
                    f"Test module: {case.cases_module}\n"
                    f"Expected: {case.unifac_result}\n"
                    f"Obtained: {result[0].subgroups}"
                )

                assert False, message


@pytest.mark.dependency(depends=["unifac"])
@pytest.mark.unifac
def test_unifac_groups_coverage():
    # Check if all the groups were detected on at least one case
    for group in TestUNIFAC.tested_groups:
        assert group in unifac.subgroups.index, f"Group {group} not tested"

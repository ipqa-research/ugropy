from cases import TCase

import pytest

from ugropy import joback


@pytest.mark.dependency(name="joback")
@pytest.mark.joback
class TestJoback(TCase):
    # Store al the groups detected in all the cases here:
    tested_groups = set()

    def asserts(self, case, solver):
        if case.joback_result is None:
            pytest.skip(
                f"No Joback result defined for {case.identifier}, "
                f"{case.cases_module}"
            )

        result = joback.get_groups(
            identifier=case.identifier,
            identifier_type=case.identifier_type,
            solver=solver,
            search_multiple_solutions=True,
        )

        if len(result) > 1:
            for r in result:
                comp = r.subgroups in case.joback_result

                if not comp:
                    message = (
                        "\n"
                        f"Case: {case.identifier}\n"
                        f"Test module: {case.cases_module}\n"
                        f"Expected: {case.joback_result}\n"
                        f"Obtained: {[r.subgroups for r in result]}"
                    )

                    assert False, message

                self.tested_groups.update(r.subgroups.keys())
        else:
            comp = result[0].subgroups == case.joback_result

            self.tested_groups.update(result[0].subgroups.keys())

            if not comp:
                message = (
                    "\n"
                    f"Case: {case.identifier}\n"
                    f"Test module: {case.cases_module}\n"
                    f"Expected: {case.joback_result}\n"
                    f"Obtained: {result[0].subgroups}"
                )

                assert False, message


@pytest.mark.dependency(depends=["joback"])
@pytest.mark.joback
def test_joback_groups_coverage():
    # Check if all the groups were detected on at least one case
    for group in TestJoback.tested_groups:
        assert group in joback.subgroups.index, f"Group {group} not tested"

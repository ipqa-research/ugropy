from cases import TCase

import pytest

from ugropy import unifac


@pytest.mark.unifac
class TestUNIFAC(TCase):
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
        else:
            comp = result[0].subgroups == case.unifac_result

            if not comp:
                message = (
                    "\n"
                    f"Case: {case.identifier}\n"
                    f"Test module: {case.cases_module}\n"
                    f"Expected: {case.unifac_result}\n"
                    f"Obtained: {result[0].subgroups}"
                )

                assert False, message

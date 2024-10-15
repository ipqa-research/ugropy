from cases import (
    alcohols_cases,
    aromathics_cases,
    hydrocarbons_cases,
    insaturated_hydrocarbons_cases,
    ketones_cases,
    particulars_cases,
    esthers_cases,
    aldehydes_cases,
)

import numpy as np


def test_no_duplicate_cases():
    # Combine all the lists into a single list
    all_cases = (
        aromathics_cases
        + alcohols_cases
        + hydrocarbons_cases
        + insaturated_hydrocarbons_cases
        + ketones_cases
        + particulars_cases
        + esthers_cases
        + aldehydes_cases
    )

    for i, case in enumerate(all_cases):
        if i == 0:
            identifiers = [case.identifier]
            continue

        if case.identifier in identifiers:
            where = np.argwhere(np.array(identifiers) == case.identifier)

            mesagge = (
                f"Duplicate found: {case.identifier}, {case.cases_module}\n"
                f"Also found in: {all_cases[where[0][0]].identifier}, "
                f"{all_cases[where[0][0]].cases_module}"
            )

            assert False, mesagge
        else:
            identifiers.append(case.identifier)

from helper import orbit, nucleon

def test_orbit():
    # Fe57
    assert nucleon('fe57')[0] == 26
    assert orbit(26) == 'pf'
    assert orbit(57-26) == 'pf'

    # Be10
    assert nucleon('be10')[0] == 4
    assert orbit(4) == '0p'
    assert orbit(10-4) == '0p'

    # O15
    assert nucleon('O15')[0] == 8
    assert orbit(8) == '0p'
    assert orbit(15-8) == '0p'

    # N15
    assert nucleon('N15')[0] == 7
    assert orbit(7) == '0p'
    assert orbit(15-7) == '0p'

    # Cu65
    assert nucleon('cu65')[0] == 29
    assert orbit(29) == 'pf'
    assert orbit(65-29) == 'pf'

    # Li7
    assert nucleon('li7')[0] == 3
    assert orbit(3) == '0p'
    assert orbit(7-3) == '0p'

    # F19
    assert nucleon('f19')[0] == 9
    assert orbit(9) == 'sd'
    assert orbit(19-9) == 'sd'

    # I127
    assert nucleon('I127')[0] == 53
    assert orbit(53) == 'jj55'
    assert orbit(127-53) == 'jj55'

    # Cs133
    assert nucleon('Cs133')[0] == 55
    assert orbit(55) == 'jj55'
    assert orbit(133-53) == 'jj55'


if __name__ == "__main__":
    try:
        test_orbit()
        print('All tests passed!!')
    except:
        print('Something wrong. Please check')

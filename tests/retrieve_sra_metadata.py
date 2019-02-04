import contextlib
import io
import os
import sys

sys.path.append(os.getcwd())

from bin.retrieve_sra_metadata import download_runs_meta



def test_missing_SRR():
    """
    Tests to make sure that if an SRR is missing it does not throw off the
    script.
    """

    with captured_output() as (out, err):
        test_set1 = ['SRR2927685', 'SRR4042625']
        download_runs_meta(test_set1)
        lines = out.getvalue().strip().split("\n")
        # There must only be one line, even though we gave to SRRs.
        assert len(lines) == 1, "test_missing_SRR: must only return 1 line."
        # We should only get the mapping for the second SRR.
        assert lines[0].strip() == "SRR4042625,SRX2033559", "test_missing_SRR: %s != %s" % ("SRR4042625,SRX2033559", lines[0].strip())



@contextlib.contextmanager
def captured_output():
    """
    Helper function to capture STDOUT/STDERR.
    """
    new_out, new_err = io.StringIO(), io.StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err



if __name__ == "__main__":
    test_missing_SRR()

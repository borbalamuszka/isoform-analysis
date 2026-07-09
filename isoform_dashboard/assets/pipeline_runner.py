import subprocess
import sys

cmd1 = ['C:\\Users\\vilsn\\Documents\\isoform-analysis\\isoform_dashboard_env\\Scripts\\python.exe', '-u', '-m', 'isoform_distribution.distributions', '--matrix', 'C:\\Users\\vilsn\\Documents\\isoform-analysis\\data\\expressed_isoforms_matrix.tsv', '--gtf', 'C:\\Users\\vilsn\\Documents\\isoform-analysis\\data\\expressed_isoforms.gtf', '--meta-file', 'C:\\Users\\vilsn\\Documents\\isoform-analysis\\data\\metadata.tsv', '--output-dir', 'data/isoform_distributions', '--meta-sample-col', 'sample_id', '--meta-group-col', 'condition', '--cutoff-pct', '1.5', '--table-type', 'both']
cmd2 = ['C:\\Users\\vilsn\\Documents\\isoform-analysis\\isoform_dashboard_env\\Scripts\\python.exe', '-u', '-m', 'isoform_distribution.bootstrap_isoform_means', '--input', 'C:\\Users\\vilsn\\Documents\\isoform-analysis\\data\\expressed_isoforms_matrix.tsv', '--output-dir', 'data/isoform_distributions', '--iterations', '1000', '--seed', '42', '--include-key', 'adult']

print("=== Running Step 1/2: Calculate Distribution Tables ===")
sys.stdout.flush()
res1 = subprocess.run(cmd1)
if res1.returncode != 0:
    print(f"\nStep 1 failed with exit code {res1.returncode}")
    sys.exit(res1.returncode)

print("\n=== Running Step 2/2: Calculate Bootstrap CIs ===")
sys.stdout.flush()
res2 = subprocess.run(cmd2)
if res2.returncode != 0:
    print(f"\nStep 2 failed with exit code {res2.returncode}")
    sys.exit(res2.returncode)

print("\n=== Preprocessing Completed Successfully! ===")

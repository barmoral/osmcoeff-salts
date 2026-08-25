jid=$(sbatch --parsable run_nacl_HP_paper-tip4p2005.sh)
echo "job 1: $jid"
for i in 2 3 4 5 6; do
    jid=$(sbatch --parsable --dependency=afterany:$jid \
          run_nacl_HP_paper-tip4p2005.sh)
    echo "job $i: $jid"
done
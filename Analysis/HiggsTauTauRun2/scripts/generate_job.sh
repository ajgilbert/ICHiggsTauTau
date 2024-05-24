if [ -z $2 ]
then
    echo "Must specify <script input> <script output>"
    exit
fi


INPUT=$1
OUTPUT=$2


echo "Generating job $OUTPUT"

echo "cd $PWD" &> $OUTPUT
echo "export X509_USER_PROXY=$HOME/cms.proxy" >> $OUTPUT
echo "source /cvmfs/cms.cern.ch/cmsset_default.sh" >> $OUTPUT
echo "export SCRAM_ARCH=slc7_amd64_gcc700" >> $OUTPUT
echo "eval \`scramv1 runtime -sh\`" >> $OUTPUT
echo "source $PWD/scripts/setup_libs.sh" >> $OUTPUT
echo "ulimit -c 0" >> $OUTPUT
echo "echo \"Cluster = \$1 Process = \$2\"" >> $OUTPUT
if [ "$INPUT" == "" ]; then :
elif [[ "$INPUT" == ./* ]];
then
  echo "$INPUT" >> $OUTPUT
else
  echo "./$INPUT" >> $OUTPUT
fi
chmod +x $OUTPUT

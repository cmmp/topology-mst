# cantor npoints seed
if [ ! -e cantor.dat ]
then
	echo "Creating data set file..."
    mvn -q exec:java -Dexec.mainClass="br.fapesp.topology.CantorSet" -Dexec.args="10000 1234" > cantor.dat
fi

# run analysis:
echo "Running analysis..."
#MAVEN_OPTS="-Xmx1024m" mvn exec:java -Dexec.mainClass="br.fapesp.topology.TopologyAnalysis" -Dexec.args="cantor.dat"
export MAVEN_OPTS="-Xmx1024m"
time mvn -q -e exec:java -Dexec.mainClass="br.fapesp.topology.TopologyAnalysis" -Dexec.args="cantor.dat" > results.dat



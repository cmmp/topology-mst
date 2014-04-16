# run analysis:
echo "Running analysis..."
#MAVEN_OPTS="-Xmx1024m" mvn exec:java -Dexec.mainClass="br.fapesp.topology.TopologyAnalysis" -Dexec.args="cantor.dat"
export MAVEN_OPTS="-Xmx1024m"
time mvn -q -e exec:java -Dexec.mainClass="br.fapesp.topology.TopologyAnalysis" -Dexec.args="logistic.dat" > results_logistic.dat



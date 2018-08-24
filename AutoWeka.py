
import subprocess

#TR == TRait
#NN == Neural Network
#RF == Random Forest
#NA == NAive
#BS == BLOSUM62
#KS == KSTAR
#SS == Random Subspace

path = "$HOME/Downloads/weka-3-8-2/weka.jar"

subprocess.run(['export CLASSPATH=' + path], shell=True)

#####SEVEN#####SEVEN#####SEVEN#####SEVEN#####SEVEN#####SEVEN#####SEVEN#####
for n in range(20):
    print(n)
    subprocess.run(["csv2arff -limit=13892 -inputs=res2,res3,res4,res6,res7,res8 bond na.csv > sevNA_balNN_" + str(n) + ".arff"], shell=True)
    subprocess.run(["java weka.classifiers.functions.MultilayerPerceptron -t sevNA_balNN_" + str(n) + ".arff -d sevNA_balNN_" + str(n) + "NN.model > sevNA_balNN_" + str(n) + "NN.out"],shell=True)


subprocess.run(["cat sevNA_balNN_*.out | grep -A 20 Stratified | grep Weighted | cat >> sevNA-NN_results.txt"],shell=True)

subprocess.run(["rm sevNA_balNN*"],shell=True)


for n in range(20):
    print(n)
    subprocess.run(["csv2arff -limit=13892 -inputs=res2,res3,res4,res6,res7,res8 bond pr.csv > sevTR_balNN_" + str(n) + ".arff"], shell=True)
    subprocess.run(["java weka.classifiers.functions.MultilayerPerceptron -t sevTR_balNN_" + str(n) + ".arff -d sevTR_balNN_" + str(n) + "NN.model > sevTR_balNN_" + str(n) + "NN.out"],shell=True)
 

subprocess.run(["cat sevTR_balNN_*.out | grep -A 20 Stratified | grep Weighted | cat >> sevTR-NN_results.txt"],shell=True)

subprocess.run(["rm sevTR_balNN*"],shell=True)


####NINE####NINE####NINE####NINE####NINE####NINE####NINE####NINE####NINE####
for n in range(20):
    print(n)
    subprocess.run(["csv2arff -limit=13892 -inputs=res1,res2,res3,res4,res6,res7,res8,res9 bond na.csv > ninNA_balNN_" + str(n) + ".arff"], shell=True)
    subprocess.run(["java weka.classifiers.functions.MultilayerPerceptron -t ninNA_balNN_" + str(n) + ".arff -d ninNA_balNN_" + str(n) + "NN.model > ninNA_balNN_" + str(n) + "NN.out"],shell=True)


subprocess.run(["cat ninNA_balNN_*.out | grep -A 20 Stratified | grep Weighted | cat >> ninNA-NN_results.txt"],shell=True)

subprocess.run(["rm ninNA_balNN*"],shell=True)



for n in range(20):
    print(n)
    subprocess.run(["csv2arff -limit=13892 -inputs=res2,res3,res4,res6,res7,res8 bond pr.csv > ninTR_balNN_" + str(n) + ".arff"], shell=True)
    subprocess.run(["java weka.classifiers.functions.MultilayerPerceptron -t ninTR_balNN_" + str(n) + ".arff -d ninTR_balNN_" + str(n) + "NN.model > ninTR_balNN_" + str(n) + "NN.out"],shell=True)
 

subprocess.run(["cat ninTR_balNN_*.out | grep -A 20 Stratified | grep Weighted | cat >> ninTR-NN_results.txt"],shell=True)

subprocess.run(["rm ninTR_balNN*"],shell=True)

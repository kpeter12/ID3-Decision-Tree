//Kevin Peterson
//OLA3
//ID3.cpp
//This program makes an ID3 decision tree based on training data
//that is provided as a command line arguement and tests that tree
//using testing data provided as a command line arguement.
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>
#include <math.h>
#define C0 1000         //these constants will be used for the attribute value of a treenode,
#define C1 1001         //they signify terminal nodes and the category which the value is classified
#define C2 1002
using namespace std;

//represents nodes in the tree
struct TreeNode
{
    int attribute;     //declares which attribute to split on, if this value is
    double value;         //the value to split on, will be left child < this value, right child > than this value
    TreeNode* left;
    TreeNode* right;
};

//function prototypes
vector<vector<int> > sort_attributes(vector<vector<double> > data);
void FindGain(vector<vector<double>> data, int numFeatures, int dataLength, TreeNode * &node);
double logBase2(double num);
vector<vector<double>> makeChildVector(double arr[][5], int dataLength);
void printTree(TreeNode * root);
void testData(double data[][5], TreeNode * root, int testDataLength, int numFeatures);

int main(int argc, char* argv[]) {

    vector<vector<double> > data;
    vector<vector<int> > indices;
    string line;
    double value;
    ifstream inFile;                        //ifstream for our training data
    ifstream inFileTesting;
    int numFeatures = atoi(argv[1]);   //number of categories of information used to classify the iris
    string trainingFile = argv[2];
    string testingFile = argv[3];
    int trainingDataLength = 0;
    int testingDataLength = 0;

    //get length of training data
    inFile.open(trainingFile);
    getline(inFile,line);
    while (!inFile.eof())
    {
        getline(inFile, line);
        trainingDataLength++;
    }
    inFile.close();
    //get length of testing data
    inFile.open(testingFile);
    getline(inFile,line);
    while (!inFile.eof())
    {
        getline(inFile, line);
        testingDataLength++;
    }
    inFile.close();
    double testArray[testingDataLength][5];    //declare array that will hold the testing data for easy manuipulation
    //cout << "Training data length: " << trainingDataLength << " Testing data length: " << testingDataLength << endl;


    inFile.open(trainingFile);

    getline(inFile,line);
    stringstream parsed(line);

    // Prep vectors...
    while (!parsed.eof()) {
        parsed >> value;
        data.push_back(vector<double>());
    }

    while (!inFile.eof()) {
        stringstream parsed(line);
        for (int i = 0; i < data.size(); i++) {
            parsed >> value;
            data[i].push_back(value);
        }
        getline(inFile,line);
    }

    //Build Tree
    TreeNode * root;
    FindGain(data,numFeatures, data[0].size(), root);
    //cout << "FINAL TREE: " << endl;
    //cout << "Root: " << root << endl;
    //printTree(root);
    inFile.close();
    //Now Test Data
    inFileTesting.open(testingFile);

    //read in testing data
    for (int i = 0; i < testingDataLength; i++) {
        for (int j = 0; j < 5; j++) {
            inFileTesting >> testArray[i][j];
        }
    }
    testData(testArray,root, testingDataLength, numFeatures);

    return 0;
}

// Attribute sorting
vector<vector<int>> sort_attributes(vector<vector<double>> data) {
    vector<vector<int> > indices;
    vector<double> *ptr;
    indices.resize(data.size());
    for (int x = 0; x < indices.size(); x++) {
        indices[x].resize(data[x].size());
        iota(indices[x].begin(),indices[x].end(),0);      //assigns all values in indices[x] to 0
        ptr = &(data[x]);
        sort(indices[x].begin(),indices[x].end(),
             [&](size_t i, size_t j){ return (*ptr)[i] < (*ptr)[j]; });
    }
    return indices;
}

//this function finds the gain and is recursively called to build the tree
void FindGain(vector<vector<double>> data, int numFeatures, int dataLength, TreeNode * &node)
{
    node = new TreeNode();
    vector<vector<double>> leftChildVector;        //vector to be given to left child, represents first half of data set
    vector<vector<double>> rightChildVector;        //vector to be given to right child, represents second half of data set
    int finalLeftChildLength, finalRightChildLength;
    int hGainAttribute = 99, hGainCutoffIndex = 99;         //highest gain attribute(index of it), highest gain cutoff
    double hGainCutoffValue = 0;                            //the value we split on (avg of hGainCutoffIndex value and the value after that)
    double hGain = -1;                        //highest gain
    double* finalArrayLeft;
    double* finalArrayRight;
    //for now split on 50, worry about multiple split points later
    vector<vector<int>> indices = sort_attributes(data);
    double dataArrayFull[dataLength][numFeatures+1];

    //first read in data to full data array and check if it is a terminal node
    for (int k = 0; k < data.size() - 1; k++) {
        for (int j = 0; j < data[0].size(); j++) {
            for (int i = 0; i < data.size(); i++) {
                dataArrayFull[j][i] = data[i][indices[k][j]];
            }
        }
    }
    //DEBUG

    /*
    cout << "Data Array FUll: " << endl;
    for (int i = 0; i < dataLength; i++) {
        for (int j = 0; j < 5; j++) {
            cout << setprecision(1) << fixed << dataArrayFull[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    */
    int numClass0 = 0, numClass1 = 0, numClass2 = 0;           //used to calc probablity of each class
    for (int i = 0; i < dataLength; i++)
    {
        if (dataArrayFull[i][numFeatures] == 0)
            numClass0++;
        else if (dataArrayFull[i][numFeatures] == 1)
            numClass1++;
        else if (dataArrayFull[i][numFeatures] == 2)
            numClass2++;
    }
    //cout << "numClass0: " << numClass0 << " numClass1: " << numClass1 << " numClass2: " << numClass2 << endl;
    if (numClass0 == dataLength)
    {
        //cout << "Node is a terminal node for C0, returning" << endl;
        node->attribute = C0;
        return;
    }
    else if (numClass1 == dataLength)
    {
        //cout << "Node is a terminal node for C1, returning" << endl;
        node->attribute = C1;
        return;
    }
    else if (numClass2 == dataLength)
    {
        //cout << "Node is a terminal node for C2, returning" << endl;
        node->attribute = C2;
        return;
    }

    //then find gain
    bool splitPoints = false;     //will check if we are able to split
    for (int k = 0; k < data.size() - 1; k++) {
        for (int count = 0; count < data[0].size(); count++) {            //for all the splits on every row.
            //will loop and try a split at each index down the row for the given attribute;
            int dataArray1RowLength = count + 1;
            int dataArray2RowLength = dataLength - (count + 1);
            double dataArray1[dataArray1RowLength][5];      //holds first half of data
            double dataArray2[dataArray2RowLength][5];      //holds second half of data

            //copies values from vector into arrays to make data easier to work with
            for (int j = 0; j < data[0].size(); j++) {
                for (int i = 0; i < data.size(); i++) {
                    if (j < count + 1) {
                        dataArray1[j][i] = data[i][indices[k][j]];
                    } else {
                        dataArray2[j - (count + 1)][i] = data[i][indices[k][j]];
                    }
                }
            }

            //only split on value changes
            if (count >= dataLength || dataArrayFull[count][k] == dataArrayFull[count + 1][k]) {
                //cout << " previous value same or index 0 or max index, skipping iteration" << endl;
                //skip iteration
            } else {
                splitPoints = true;
                //FIND GAIN HERE, Current attribute = k,
                //E(attribute) = P(S1) <- which is the dataArray1 length
                double pS1 = double(dataArray1RowLength) / double(dataLength);
                double pS2 = double(dataArray2RowLength) / double(dataLength);
                //cout << "pS1: " << pS1 << " pS2: " << pS2 << endl;

                //find the probability of all the classes given S1
                double class0CountGivenS1 = 0, class1CountGivenS1 = 0, class2CountGivenS1 = 0;
                for (int y = 0; y < dataArray1RowLength; y++) {
                    //cout << "dataArray[" <<y <<"][numFeatures]: "<<dataArray1[y][numFeatures] << "     " << dataArray1[y][4] << endl;
                    if (dataArray1[y][numFeatures] == 0)
                        class0CountGivenS1++;
                    else if (dataArray1[y][numFeatures] == 1)
                        class1CountGivenS1++;
                    else if (dataArray1[y][numFeatures] == 2)
                        class2CountGivenS1++;
                }
                class0CountGivenS1 /= dataArray1RowLength;
                class1CountGivenS1 /= dataArray1RowLength;
                class2CountGivenS1 /= dataArray1RowLength;
                //find the probability of all the classes given S2
                double class0CountGivenS2 = 0, class1CountGivenS2 = 0, class2CountGivenS2 = 0;
                for (int y = 0; y < dataArray2RowLength; y++) {
                    if (dataArray2[y][numFeatures] == 0)
                        class0CountGivenS2++;
                    else if (dataArray2[y][numFeatures] == 1)
                        class1CountGivenS2++;
                    else if (dataArray2[y][numFeatures] == 2)
                        class2CountGivenS2++;
                }
                class0CountGivenS2 /= dataArray2RowLength;
                class1CountGivenS2 /= dataArray2RowLength;
                class2CountGivenS2 /= dataArray2RowLength;

                double expectedInfo = double(-1) * pS1 * (class0CountGivenS1 * logBase2(class0CountGivenS1) +
                                                          class1CountGivenS1 * logBase2(class1CountGivenS1)
                                                          + class2CountGivenS1 * logBase2(class2CountGivenS1)) -
                                      pS2 * (class0CountGivenS2 * logBase2(class0CountGivenS2) +
                                             class1CountGivenS2 * logBase2(class1CountGivenS2) +
                                             class2CountGivenS2 * logBase2(class2CountGivenS2));
                //cout << "E(" << k << ") = " << setprecision(5) << expectedInfo << endl;

                //I(0,1,2) = H(0) + H(1) + H(2)
                //first find probabilities
                double probClass0 = 0, probClass1 = 0, probClass2 = 0;
                for (int z = 0; z < dataLength; z++) {
                    if (dataArrayFull[z][numFeatures] == 0)
                        probClass0++;
                    else if (dataArrayFull[z][numFeatures] == 1)
                        probClass1++;
                    else if (dataArrayFull[z][numFeatures] == 2)
                        probClass2++;
                }
                probClass0 /= dataLength;
                probClass1 /= dataLength;
                probClass2 /= dataLength;
                //now plug into equation
                double totalInformation =
                        -1 * probClass0 * (logBase2(probClass0)) - probClass1 * (logBase2(probClass1)) -
                        probClass2 * (logBase2(probClass2));
                //cout << "I(0,1,2) = " << setprecision(5) << totalInformation << endl;

                //gain = I(0,1,2) - E(K)
                double gain = totalInformation - expectedInfo;
                //cout << "Gain = " << setprecision(5) << gain << endl;

                //update gain
                if ((gain > hGain) || (gain == hGain && k < hGainAttribute) || ((gain == hGain) && (k <= hGainAttribute) &&
                                                                                       ((dataArrayFull[hGainCutoffIndex][hGainAttribute] +
                                                                                       dataArrayFull[hGainCutoffIndex + 1][hGainAttribute]) / 2) < hGainCutoffValue)) {
                    hGain = gain;
                    hGainAttribute = k;
                    hGainCutoffIndex = count;  //cutoff after this value
                    hGainCutoffValue = (dataArrayFull[hGainCutoffIndex][hGainAttribute] +
                                        dataArrayFull[hGainCutoffIndex + 1][hGainAttribute]) / 2;
                    //cout << "Gain updated...  hGain=" << hGain << " hGainAttribute=" << k << " hGainCutoffIndex="
                    //     << hGainCutoffIndex << " hGainCutoffValue = " << hGainCutoffValue << endl;

                    leftChildVector = makeChildVector(dataArray1,dataArray1RowLength);
                    rightChildVector = makeChildVector(dataArray2,dataArray2RowLength);
                    finalLeftChildLength = dataArray1RowLength;
                    finalRightChildLength = dataArray2RowLength;

                    node->attribute = hGainAttribute;
                    node->value = hGainCutoffValue;

                    //cout << " node attribute: " << node->attribute << " node value: " << node->value << endl;
                } else {
                    //cout << "Gain not updated" << endl;
                }

                //END FIND GAIN
                /*
                cout << "printing dataArray1 for count" << count << ": " << endl;
                for (int i = 0; i < count + 1; i++) {
                    for (int j = 0; j < 5; j++) {
                        cout << setprecision(1) << fixed << dataArray1[i][j] << " ";
                    }
                    cout << endl;
                }

                cout << "printing dataArray2 for count" << count << ": " << endl;
                for (int i = 0; i < 100 - (count + 1); i++) {
                    for (int j = 0; j < 5; j++) {
                        cout << setprecision(1) << fixed << dataArray2[i][j] << " ";
                    }
                    cout << endl;
                }
                */
                //cout << endl << "Final hGain=" << hGain << " hGainAttribute=" << k << " hGainCutoffIndex="
                //             << hGainCutoffIndex << " hGainCutoffValue = " << hGainCutoffValue << endl;

            }
        }

    }
    //Check for ties
    if (!splitPoints)
    {
        //cout << "TIE" << endl;
        //cout << "DEBUG: Making extra split" << endl;
        int count0 = 0, count1 = 1, count2 = 0;
        for (int i = 0; i < dataLength; i++)
        {
            if (dataArrayFull[i][numFeatures] == 0)
                count0++;
            else if (dataArrayFull[i][numFeatures] == 1)
                count1++;
            else if (dataArrayFull[i][numFeatures] == 2)
                count2++;
        }
        if (count0 >= count1) {
            node->attribute = C0;
            return;
        }
        else if (count1 >= count2){
            node->attribute = C1;
            return;
        }
        else{
            node->attribute = C2;
            return;
        }
    }
    //cout << "Final Node attribute: " << node->attribute << " final node value: " << node->value << endl;
    FindGain(leftChildVector,numFeatures,finalLeftChildLength,node->left);
    FindGain(rightChildVector,numFeatures,finalRightChildLength,node->right);
    //cout << "root: " << node << endl;

}

//works like regualar log base 2 function but returns 0 where log base two returns nan
double logBase2(double num)
{
    if (num == 0)
    {
        //cout << " num was zero, returning 0 for log  ";
        return 0;
    }

    else if (num < 0)
    {
        //cout << "log base 2 num: " << num << " is below 0" << endl;
        return 0;
    }
    else
        return log2(num);
}
//makes vectory out of array to pass down to recursive child function
vector<vector<double>> makeChildVector(double arr[][5], int dataLength)
{
    vector<vector<double>> childVector;
    double value;
    // Prep vectors...
   for(int i = 0; i < 5; i++) {
        value = arr[0][i];
        childVector.push_back(vector<double>());
    }

    for (int i = 0; i < dataLength; i++) {
        for (int j = 0; j < childVector.size(); j++) {
            value = arr[i][j];
            childVector[j].push_back(value);
        }
    }

    return childVector;
}
//recursively prints tree for debugging
void printTree(TreeNode * root)
{
    if (root == NULL)
    {
        return;
    }
    else
    {
        cout << root->attribute << " ";
        cout << " < "<<root->value <<" ";
        printTree(root->left);
        cout << " > " << root->value << " ";
        printTree(root->right);
    }
}

//will print the number of data points that the tree correctly classified
void testData(double data[][5], TreeNode * root, int testDataLength, int numFeatures)
{
    int numCorrect = 0;
    for (int i = 0; i < testDataLength; i++)
    {
        TreeNode * node = root;
        //cout << "testinng line: " << data[i][0] << " " << data[i][1] << " "<< data[i][2] << " "<< data[i][3] << " "<< data[i][4] << " " << endl;
        int correctCategory = int(data[i][numFeatures] + 1000);
        while (node->attribute < 999) //while not a terminal node
        {
            if (data[i][node->attribute] <= node->value)
            {
                //cout << "     data: " << data[i][node->attribute] << "  value: " << node->value << "going to left child" << endl;
                node = node->left;
            }
            else if (data[i][node->attribute] > node->value)
            {
                //cout << "     data: " << data[i][node->attribute] << "  value: " << node->value << "going to right child" << endl;
                node = node->right;
            }
        }
        if (correctCategory == node->attribute)
        {
            //cout << "Tree correctly classified this value   value: " << node->attribute << "  correct category: "<< correctCategory << endl;
            numCorrect++;
        }
        else
        {
            //cout << "Tree classified wrong value" << node->attribute << "  correct category: "<< correctCategory << endl;
        }
    }
    cout << numCorrect << endl;
}

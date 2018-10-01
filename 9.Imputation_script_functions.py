#=============================================================================#
#================================ LIBRARIES ==================================#
#=============================================================================#

print("\n\n#################\nLoading Libraries")
## MAIN ##
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## PIPELINE ##
from sklearn.pipeline import Pipeline
from sklearn.cross_validation import StratifiedKFold
from sklearn.grid_search import GridSearchCV

# from sklearn import  metrics
from sklearn.metrics import accuracy_score

import sklearn.metrics 
from sklearn.preprocessing import StandardScaler, MinMaxScaler

## FEATURE ##
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import Lasso
from sklearn.feature_selection import VarianceThreshold

## CLASSIFIERS ##
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.feature_selection import SelectPercentile,f_classif
transform = SelectPercentile(f_classif)
from sklearn.model_selection import cross_val_score
from sklearn.metrics import confusion_matrix
from sklearn.decomposition import PCA
from sklearn.multiclass import OneVsRestClassifier
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from scipy import interp



#=====================================================================#
#============================== LD DICTIONARY ========================#
def ADD_TO_LD_DICT(name1,name2,LD_dict):
    if name1 in LD_dict.keys():
        LD_dict[name1].append(name2)
    else:
        LD_dict[name1] = [name2]
    if name2 in LD_dict.keys():
        LD_dict[name2].append(name1)
    else:
        LD_dict[name2] = [name1]
            
    return LD_dict
#--------------------------------------------------------------------#

#=============================================================================#
#============================== Choose Classifier ============================#
#=============================================================================#
#clas = input("Choose Classifier\n1. KNN\n2.Logistic Regression\n3. SVM\n((1/2/3) = ")
# print(clas+"ok")
def Classifier_Ch(choice_cl,Min_F,thres,S,D,num_clas):
    '''
    **Description:**\n 
    Create a pipeline with stdscale, variance feature selection and classifier (user input).\n
    **Input:**\n
    - Classifier: 1 = KNN, 2 = Logistic Regression, 3. SVM\n
    - Min_f: Number of Folds\n
    - Variance threshold\n
    **Output:**\n
    - Pipeline 
    '''
    ## Assign Feature selection ##
    # selector = VarianceThreshold(threshold = 0.001)
    ## Assign Classifier ##
    if choice_cl == "KNN" :
        clf = KNeighborsClassifier()
        param_grid = [{'clf__n_neighbors': list(range(1, 20,2)),'clf__metric': ['minkowski','euclidean','manhattan'] ,'clf__weights':['uniform','distance']}]
    elif choice_cl == "LG" :
        clf = LogisticRegression()
        param_grid = [{'clf__penalty': ['l1','l2'],'clf__C': [0.1, 10, 100]}]
    elif choice_cl == "RF" :
        clf = RandomForestClassifier()
        param_grid = [{'clf__max_features': ['auto',None],'clf__bootstrap': [True],'clf__max_depth':[None,2,4,5,6]}]
    else:
        clf = SVC()
        param_grid = [{'clf__kernel': ['linear'],'clf__C': [0.01, 0.1, 10, 100],'clf__probability':[True]},{'clf__kernel': ['poly'],'clf__C': [0.0001, 0.001, 0.01, 0.1, 10, 100],'clf__degree' :[2,3],'clf__probability':[True]}]
        if num_clas==2:
            param_grid = [{'clf__kernel': ['linear'],'clf__C': [0.01, 0.1, 10, 100],'clf__probability':[True],'clf__decision_function_shape':['ovr']},{'clf__kernel': ['poly'],'clf__C': [0.0001, 0.001, 0.01, 0.1, 10, 100],'clf__degree' :[2,3],'clf__probability':[True],'clf__decision_function_shape':['ovr']}]

    ## Create Pipeline ##
    pipe = Pipeline([('clf', clf)])
    gcv = GridSearchCV(estimator=pipe,param_grid=param_grid,scoring='accuracy',cv=Min_F)
    # gcv = GridSearchCV(estimator=pipe,param_grid=param_grid,scoring='roc_auc',cv=Min_F)
    return gcv,choice_cl




#=========================================================================
#========================= CROSS VALIDATION ==============================
#=========================================================================
def CV_Kyriakis(SNP,Train_Data,Known_Labels,New_Data,name,Min_Folds,Test_Data_samples):
    
    uni_clas  =len(np.unique(Known_Labels)) 
    

    thres = 0.01
    in_out = "IN"

    mat = np.zeros((2,2))
    mean_acc = 0
    Best_Param = {}
    lista_acc =[]
    ### REPEATED NESTED ###
    # output = open("Results.txt","w")
    for  i in [1]:
        ## NESTED CROSS VALIDATION ##
        kfold = StratifiedKFold(y=Known_Labels, n_folds= Min_Folds, shuffle=True)
       
        svm_list = []
        knn_list =[]
        lg_list = []
        
        for train_idx, test_idx in kfold:
            ## Assign ##
            inner_Train = Train_Data[train_idx]
            S = inner_Train.shape[0]
            D = inner_Train.shape[1]
            inner_Test = Train_Data[test_idx]
            inner_labels = Known_Labels[train_idx]

            #============ CHECK NUM FOLDS ============#
            
            Num_HomR = sum(inner_labels == "0")
            Num_Hetr = sum(inner_labels == "1")
            Num_HomA = sum(inner_labels == "2")

            
            # Set minimum Folds Based on classes
            Min_Folds_inner = min(Num_HomR,Num_Hetr,Num_HomA)
            if Num_HomA == 0:
                Min_Folds_inner = min(Num_HomR,Num_Hetr)
            if Num_HomR == 0:
                Min_Folds_inner = min(Num_Hetr,Num_HomA)
            if Num_Hetr == 0:
                Min_Folds_inner = min(Num_HomR,Num_HomA)
            #---------------------------------------------------------------------#
            if Min_Folds_inner ==2 or Min_Folds_inner ==1:
                continue
            if Min_Folds_inner >=5:
                Min_Folds_inner = 5
            
            ## Prepare ##
            gcv,name = Classifier_Ch(name ,Min_Folds_inner,thres,S,D,uni_clas)     
            ## CROSS VALIDATION ##
            gcv.fit(inner_Train,inner_labels)
            best_par = gcv.best_params_
           
            if SNP in Best_Param.keys():
                Best_Param[str(best_par)]['counter'] +=1
            else:
                best_par['counter'] = 1
                Best_Param[str(best_par)] = best_par
            ## Predict ##
            y_test = Known_Labels[test_idx]
            y_pred = gcv.predict(inner_Test)
            
            ## Accuracy ##
            acc = accuracy_score(y_true=y_test, y_pred=y_pred)
            lista_acc.append(acc)
            # output.write("\n"+ str(best_par)+' | inner ACC %.2f%% | outer ACC %.2f%%' % (gcv.best_score_ * 100, acc * 100))
            mean_acc = mean_acc + acc

    mean_acc_f =mean_acc/(1*Min_Folds)
    # output.write(str(mean_acc_f))
    # output.write("\n"+str(lista_acc))
    # for x in Best_Param:
    #     output.write("\n"+str(x)+" : "+str(Best_Param[x]))
    # output.close()
    
    if mean_acc_f==0.0:
        pipe=0
        predict=0
        return mean_acc_f,pipe,predict
    # print(str(mean_acc_f))
    Best = max(Best_Param, key=lambda i: Best_Param[i]['counter'])
    
    # print(Best)
    if name=='LG':
        C = Best_Param[Best]['clf__C']
        PEN = Best_Param[Best]['clf__penalty']
        pipe = Pipeline([('clf', LogisticRegression(C=C,penalty=PEN))])
    elif name=='SVM':
        C = Best_Param[Best]['clf__C']
        kernel = Best_Param[Best]['clf__kernel']

        pipe = Pipeline([('clf', SVC(C=C,kernel=kernel, probability=True))])
        # ADD MULTICLASS OVR
        if uni_clas ==2:
            pipe = Pipeline([('clf', SVC(C=C,kernel=kernel, probability=True,decision_function_shape='ovr'))])
        
        if 'clf__degree' in Best_Param[Best].keys():
            degree = Best_Param[Best]['clf__degree']
            pipe = Pipeline([('clf', SVC(C=C,kernel=kernel,degree=degree, probability=True))])
            if uni_clas==2:
                pipe = Pipeline([('clf', SVC(C=C,kernel=kernel,degree=degree, probability=True,decision_function_shape='ovr'))])
    elif name=='KNN':
        neighbors = Best_Param[Best]['clf__n_neighbors']
        metric = Best_Param[Best]['clf__metric']
        weights = Best_Param[Best]['clf__weights']
        pipe = Pipeline([('clf', KNeighborsClassifier(metric = metric, n_neighbors = neighbors, weights = weights))])
    elif name=='RF':
        max_features = Best_Param[Best]['clf__max_features']
        max_depth = Best_Param[Best]['clf__max_depth']
        pipe = Pipeline([('clf', RandomForestClassifier(bootstrap=True,max_features=max_features,max_depth=max_depth))])

    #=================================== Predict =================================#
    pipe.fit(Train_Data,Known_Labels)
    ### 4.  Predict ##
    predictions  = pipe.predict(New_Data)
    ## 5.  Write ##
    predict = {'Id': Test_Data_samples, 'Predicted':list(predictions)}

    return mean_acc_f,pipe,predict    # print(predict)
    # predict1 = pd.DataFrame(predict)
    # predict1.to_csv("Predicted.csv", index = False)
#------------------------------------------------------------------------#

#=========================================================================
#================================ ROC - AUC ==============================
#=========================================================================

def ROC_Kyriakis(SNP,Train_Data,Known_Labels,name,mean_acc_f,predict,ld,Min_Folds,pipe,pos_label):
   
    from sklearn.linear_model import Lasso
    cv = StratifiedKFold(Known_Labels, n_folds=Min_Folds,random_state=1)
    fig = plt.figure(figsize=(10, 8))
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_tpr = []
    

    for i, (train, test) in enumerate(cv):
        probas = pipe.fit(Train_Data[train],Known_Labels[train]).predict_proba(Train_Data[test])
        fpr, tpr, thresholds = roc_curve(Known_Labels[test], probas[:, 1],pos_label=pos_label)
        if pos_label==0:
            tpr = 1-tpr
            fpr = np.short(1-fpr)
        mean_tpr += interp(mean_fpr, fpr, tpr)

        mean_tpr[0] = 0.0
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr,tpr,lw=1,label='ROC fold %d (area = %0.2f)'% (i+1, roc_auc))
    plt.plot([0, 1],[0, 1],linestyle='--',color=(0.6, 0.6, 0.6),label='random guessing')
    mean_tpr /= len(cv)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    if mean_auc>=0.85:
        plt.plot(mean_fpr, mean_tpr, 'k--', label='mean ROC (area = %0.2f)' % mean_auc, lw=2)
        plt.plot([0, 0, 1], [0, 1, 1],lw=2,linestyle=':',color='black',label='perfect performance')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('false positive rate')
        plt.ylabel('true positive rate')
        plt.title('Receiver Operator Characteristic\n'+SNP+"\n"+str(predict)+"\nAccuracy = {0:.3f}".format(mean_acc_f))
        plt.legend(loc="lower right")
        plt.savefig("Impute_Images/"+SNP.replace(":","_")+'_'+name+'_{}_{}.png'.format(ld,pos_label))
    plt.close()
    return mean_auc




def ROC_VS(SNP,Train_Data,Known_Labels,name,mean_acc_f,predict,ld,Min_Folds,pipe):
    from sklearn.preprocessing import label_binarize
    from sklearn.model_selection import train_test_split
    from itertools import cycle
    n_classes = len(np.unique(Known_Labels))
    # print(n_classes)
    # clas = np.unique(Known_Labels).tolist()
    # print(clas)
    # Add noisy features to make the problem harder
    random_state = np.random.RandomState(0)
    
    # Learn to predict each class against the other
    classifier = OneVsRestClassifier(pipe)
    # svm.SVC(kernel='linear', probability=True,random_state=random_state))
    
    cv = StratifiedKFold(Known_Labels.T, n_folds=Min_Folds,random_state=0)
    micro = 0
    macro = 0

    prev_macro = 0 
    for ji, (train, test) in enumerate(cv):
        X_train = Train_Data[train]
        X_test  = Train_Data[test]
        
        y = label_binarize(Known_Labels, classes=[0,1,2])
        # nums = y.shape[1]
        
        y_train = y[train,:]
        y_test  = y[test,:]

        y_score = classifier.fit(X_train, y_train).decision_function(X_test)
        # Compute ROC curve and ROC area for each class
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

        # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])    

        # =================== PLOT ========================== #
        # Compute macro-average ROC curve and ROC area

        # First aggregate all false positive rates
        all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

        # Then interpolate all ROC curves at this points
        mean_tpr = np.zeros_like(all_fpr)
        for i in range(n_classes):
            mean_tpr += interp(all_fpr, fpr[i], tpr[i])

        # Finally average it and compute AUC
        mean_tpr /= n_classes

        fpr["macro"] = all_fpr
        tpr["macro"] = mean_tpr
        roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])
        
        micro += roc_auc["micro"]
        macro += roc_auc["macro"]
        
        # Plot all ROC curves
        if  roc_auc["macro"]>prev_macro:
            fpr_Mi, tpr_Mi = fpr["micro"], tpr["micro"]
            fpr_Ma, tpr_Ma = fpr["micro"], tpr["micro"]
            fpr_f , tpr_f  = fpr,tpr
            max_macro = roc_auc["macro"]
            max_micro = roc_auc["micro"]
            roc_auc_f = roc_auc
        prev_macro =  roc_auc["macro"]
    mean_micro = micro/Min_Folds
    mean_macro = macro/Min_Folds
    if mean_micro>=0.85 and mean_macro>=0.85:

        fig = plt.figure(figsize=(10, 8))
        plt.plot(fpr_Mi, tpr_Mi,
                 label='micro-average ROC curve (area = {0:0.2f})'
                       ''.format(max_micro),
                 color='deeppink', linestyle=':', linewidth=4)

        plt.plot(fpr_Ma, tpr_Ma,
                 label='macro-average ROC curve (area = {0:0.2f})'
                       ''.format(max_macro),
                 color='navy', linestyle=':', linewidth=4)

        colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
        for i, color in zip(range(n_classes), colors):
            plt.plot(fpr_f[i], tpr_f[i], color=color, lw=2,
                     label='ROC curve of class {0} (area = {1:0.2f})'
                     ''.format(i, roc_auc_f[i]))

        plt.plot([0, 1], [0, 1], 'k--', lw=2)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver operating characteristic to multi-class\n{}\n{} fold Micro,Macro = {:.3f},{:.3f}\n{}'.format(SNP,Min_Folds,mean_micro,mean_macro,str(predict)))
        plt.legend(loc="lower right")
        plt.savefig("Impute_Images/"+SNP.replace(":","_")+'_'+name+'_{}.png'.format(ld))
        plt.close()
    return(mean_micro,mean_macro)




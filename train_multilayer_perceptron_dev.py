# ----------------------------------------------------------------------------
# train_multilayer_perceptron.py 
# Author: Antoine Beauchamp
# Created: December 15th, 2020

"""
Train a multi-layer perceptron

Description
-----------

"""

# Packages -------------------------------------------------------------------

import argparse
import os
import random
import pandas as pd
import numpy as np
from datatable import fread

import torch
import torch.nn.functional as F
from torch                    import nn
from torch.optim              import AdamW
from torch.optim.lr_scheduler import OneCycleLR
from torch.cuda               import is_available

from skorch                   import NeuralNetClassifier
from skorch.callbacks         import LRScheduler
from skorch.helper            import DataFrameTransformer

from sklearn.metrics import accuracy_score, confusion_matrix

# Functions ------------------------------------------------------------------

def parse_args():
    
    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )
    
    parser.add_argument(
        "--datadir",
        type = str,
        default = 'data/',
        help = "Directory containing input data."
    )
    
    parser.add_argument(
        "--outdir",
        type = str,
        default = 'data/MLP_outcomes/',
        help = "Directory in which to write neural net outcomes."
    )
    
    parser.add_argument(
        '--training',
        type = str,
        help = "Path to CSV file containing training data."
    )
    
    parser.add_argument(
        '--mousetransform',
        type = str,
        help = ("Path to CSV file containing mouse data to transform")
    )
    
    parser.add_argument(
        '--humantransform',
        type = str,
        help = ("Path to CSV file containing human data to transform")
    )
    
#     parser.add_argument(
#         "--labels",
#         type = str,
#         default = 'region5',
#         choices = ['region5', 
#                    'region11',
#                    'region28',
#                    'region46',
#                    'region67',
#                    'region134'],
#         help = "Class of mouse labels on which to train."
#     )
    
#     parser.add_argument(
#         "--mousedata",
#         type = str,
#         default = 'region134',
#         choices = ['region5', 
#                    'region11',
#                    'region28',
#                    'region46',
#                    'region67',
#                    'region134'],
#         help = "Mouse data to apply the trained network to."
#     )
    
#     parser.add_argument(
#         "--humandata",
#         type = str,
#         default = 'region166',
#         choices = ['region5',
#                    'region16',
#                    'region56',
#                    'region79',
#                    'region88',
#                    'region166'],
#         help = "Human data to apply the trained network to."
#     )
    
    parser.add_argument(
        "--nunits",
        type = int,
        default = 500,
        help = "Number of hidden units in the network."
    )
    
    parser.add_argument(
        "--L2",
        type = float,
        default = 1e-6,
        help = "Weight decay"
    )
    
    parser.add_argument(
        "--nepochs",
        type = int,
        default = 200,
        help = "Number of epochs to train over."
    )
    
    parser.add_argument(
        "--learningrate",
        type = float,
        default = 1e-5,
        help = "Learning rate during training."
    )
    
    parser.add_argument(
        "--confusionmatrix",
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = ("Option to compute confusion matrix from training set "
                "predictions.")
    )
    
    parser.add_argument(
        "--voxeltransform",
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = ("Option to transform the voxel-wise data using the "
                "modified MLP.")
    )
    
    parser.add_argument(
        '--seed',
        type = int,
        help = ("Random seed")
    )
    
    parser.add_argument(
        '--verbose',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = ("Verbosity.")
    )
        
    
    args = vars(parser.parse_args())
    
    return args

def prepare_data(data, labelcol):

    """
    """

    ind_labels = data.columns.str.match(labelcol)

    input_data = data.loc[:, ~ind_labels].copy()
    label_data = data.loc[:, ind_labels].copy()

    label_data.loc[:, labelcol] = label_data.loc[:, labelcol].astype('category')

    X = DataFrameTransformer().fit_transform(input_data)['X']
    y = DataFrameTransformer().fit_transform(label_data)[labelcol]

    return X, y


def define_classifier(input_units, output_units, hidden_units, weight_decay, max_epochs, learning_rate, train_split = None):

    #Define network architecture
    class ClassifierModule(nn.Module):
        def __init__(
            self,
            input_units,
            output_units,
            hidden_units,
            apply_output_layer = True #Flag to apply output layer
        ):
            super(ClassifierModule, self).__init__()

            self.apply_output_layer = apply_output_layer

            self.hidden1 = nn.Linear(input_units, hidden_units)
            self.hidden2 = nn.Linear(hidden_units, hidden_units)
            self.hidden3 = nn.Linear(hidden_units, hidden_units)
            self.output = nn.Linear(hidden_units, output_units)

        def forward(self, X, **kwargs):
            X = F.relu(self.hidden1(X))
            X = F.relu(self.hidden2(X))
            X = F.relu(self.hidden3(X))

            #If flag is True, apply output layer
            if self.apply_output_layer is True:
                X = F.softmax(self.output(X), dim = -1)

            return X

    #Create the classifier
    net = NeuralNetClassifier(
            ClassifierModule(input_units = input_units,
                             output_units = output_units,
                             hidden_units = hidden_units),
            train_split = None,
            optimizer = AdamW,
            optimizer__weight_decay = weight_decay,
            max_epochs = max_epochs,
            callbacks = [('lr_scheduler',
                          LRScheduler(policy=OneCycleLR,
                                      total_steps=max_epochs,
                                      cycle_momentum=False,  
                                      max_lr=learning_rate))] 
        )

    return net


def calc_confusion_matrix(data, y_true, y_pred, labelcol):

    """
    """

    df_confmat = pd.DataFrame(confusion_matrix(y_true, y_pred))

    df_labels = data[[labelcol]].copy()
    df_labels['Dummy'] = y_true
    df_labels = df_labels.sort_values('Dummy').drop_duplicates().reset_index(drop = True)

    df_confmat.columns = df_labels[labelcol].astype('str')

    df_confmat['TrueLabels'] = (df_labels[labelcol]
                                .astype('str')
                                .reset_index(drop = True))

    return df_confmat
    
    
# Main ------------------------------------------------------------------------

def main():

    #Load command line arguments
    args = parse_args()
    
    datadir = args['datadir']
    outdir = args['outdir']
    confmat = True if args['confusionmatrix'] == 'true' else False
    verbose = True if args['verbose'] == 'true' else False
    
    datadir = os.path.join(datadir, '')
    outdir = os.path.join(outdir, '')
    
    if os.path.exists(outdir) == False:
        print('Output directory {} not found. Creating it...'.format(outdir))
        os.mkdir(outdir)
    
    
    # Import data -------------------------------------------------------------

    #Set up files for import
    #Mouse voxelwise data to train over
    training_file = args['training']
    
    if verbose:
        print("Importing data...")

    #Import data
    df_training = (fread(training_file, header = True)
                   .to_pandas())

    
    # Process data ------------------------------------------------------------

    if verbose:
        print("Preparing data for learning...")
    
    X, y = prepare_data(data = df_training, 
                        labelcol = 'Region')
    

    # Train the network ------------------------------------------------------

    if verbose:
        print("Initializing neural network...")

    seed = args['seed']
    if seed is not None:
        np.random.seed(seed)
        torch.manual_seed(seed)
        random.seed(seed)
        
    #Get network parameters from command line args
    hidden_units = args['nunits']
    weight_decay = args['L2']
    max_epochs = args['nepochs']
    learning_rate = args['learningrate']
    nlabels = len(np.unique(y))
    
    net = define_classifier(input_units = X.shape[1],
                            output_units = nlabels,
                            hidden_units = hidden_units,
                            weight_decay = weight_decay,
                            max_epochs = max_epochs,
                            learning_rate = learning_rate)
    
    if is_available() == True:
        print("GPU available. Training network using GPU...")
    else:
        print("GPU unavailable. Training network using CPU...")
    
    #Fit the network
    net.fit(X, y)
    
    #Predict training labels
    y_pred = net.predict(X)
    
    #Compute training accuracy
    if verbose:
        print("Training accuracy: {}".format(accuracy_score(y, y_pred)))
    
    
    # Compute training confusion matrix --------------------------------------
 
    if confmat:
        
        if verbose:
            print("Computing confusion matrix from training set...")
        
        df_confmat = calc_confusion_matrix(data = df_training,
                                           y_true = y,
                                           y_pred = y_pred,
                                           labelcol = 'Region')
        
        #File to save confusion matrix
        file_confmat = "MLP_confusionmatrix_training"+\
        "_labels"+str(nlabels)+\
        "_layers3"+\
        "_units"+str(args['nunits'])+\
        "_L2"+str(args['L2'])+".csv"
        
        #Write confusion matrix to file
        df_confmat.to_csv(os.path.join(outdir, file_confmat), 
                          index = False)


    # Predict label probabilities for mouse/human data -----------------------
    
    def transform_input_space(data, network, labelcol):
        
        X, y = prepare_data(data = data, 
                            labelcol = labelcol)
        
        df_transformed = pd.DataFrame(network.predict_proba(X))
        
        df_transformed[labelcol] = data[labelcol]
        
        return df_transformed
    
    mouse_transform = args['mousetransform']
    human_transform = args['humantransform']
    
    print(mouse_transform)
    print(human_transform)
    
    if mouse_transform is not None:
        df_mouse = (fread(mouse_transform, header = True)
                    .to_pandas())
    else: 
        df_mouse = df_training


    df_mouse_prob = transform_input_space(data = df_mouse,
                                          network = net,
                                          labelcol = 'Region')
    
    df_mouse_prob.columns = np.append(np.unique(df_mouse['Region']), 'Region')
    
    nlabels_mouse = len(np.unique(df_mouse['Region']))
    
    #File to save mouse probabilities
    file_mouse_prob = 'MLP'+\
    '_labels'+str(nlabels)+\
    '_layers3'+\
    '_units'+str(args['nunits'])+\
    '_L2'+str(args['L2'])+\
    '_mouseprob_'+str(nlabels_mouse)+\
    '.csv'
    
    print(file_mouse_prob)
    quit()
    
    # STOPPED HERE ---------------
    
    
    if human_transform is not None:
        
        df_human = (fread(human_transform, header = True)
                    .to_pandas())
        
        df_human_prob = transform_input_space(data = df_human,
                                              network = net,
                                              labelcol = 'Region')
        
        df_human_prob.columns = np.append(np.unique(df_human['Region']), 'Region')
        
    
    
    quit()
        
    
    if (mouse_transform is None) & (human_transform is None):
        print("Test")
    
        
    print("Applying trained network to mouse and human data...")

        
    
    #Put mouse/human data into appropriate format for network        
    X_Mouse = dftx.fit_transform(dfInputMouse)['X']
    X_Human = dftx.fit_transform(dfInputHuman)['X']
    
    
    #Compute label probabilities for mouse/human data
    dfPredictionsMouse = pd.DataFrame(net.predict_proba(X_Mouse))
    dfPredictionsHuman = pd.DataFrame(net.predict_proba(X_Human))
    
    #Include label names as columns
    dfPredictionsMouse.columns = dfLabelsUnique[labelcol].astype('str')
    dfPredictionsHuman.columns = dfLabelsUnique[labelcol].astype('str')
    
    #Include true labels
    dfPredictionsMouse['TrueLabel'] = dfExprMouse['Region']
    dfPredictionsHuman['TrueLabel'] = dfExprHuman['Region']
    
    #File to save mouse probabilities
    fileMouseProb = "MLP_"+args['labels'].capitalize()+\
    "_Layers3"+\
    "_Units"+str(args['nunits'])+\
    "_L2"+str(args['L2'])+\
    "_MouseProb_"+args['mousedata'].capitalize()+'.csv'
    
    #File to save human probabilities
    fileHumanProb = "MLP_"+args['labels'].capitalize()+\
    "_Layers3"+\
    "_Units"+str(args['nunits'])+\
    "_L2"+str(args['L2'])+\
    "_HumanProb_"+args['humandata'].capitalize()+'.csv'
    
    #Write probability data to file
    dfPredictionsMouse.to_csv(os.path.join(outdir, fileMouseProb), 
                              index = False)
    dfPredictionsHuman.to_csv(os.path.join(outdir, fileHumanProb), 
                              index = False)
    
    # Extract hidden layer for mouse/human data ------------------------------

    #Change the mode of the network to extract a hidden layer
    net.module_.apply_output_layer = False
    
    #Apply the modified network to the mouse and human data to get the
    #hidden units
    dfMouseTransformed = pd.DataFrame(net.predict_proba(X_Mouse))
    dfHumanTransformed = pd.DataFrame(net.predict_proba(X_Human))
    
    #Include region information
    dfMouseTransformed['Region'] = dfExprMouse['Region']
    dfHumanTransformed['Region'] = dfExprHuman['Region']

    #File to save transformed mouse data
    fileMouseTx = "MLP_"+args['labels'].capitalize()+\
    "_Layers3"+\
    "_Units"+str(args['nunits'])+\
    "_L2"+str(args['L2'])+\
    "_MouseTx_"+args['mousedata'].capitalize()+'.csv'
    
    #File to save transformed human data
    fileHumanTx = "MLP_"+args['labels'].capitalize()+\
    "_Layers3"+\
    "_Units"+str(args['nunits'])+\
    "_L2"+str(args['L2'])+\
    "_HumanTx_"+args['humandata'].capitalize()+'.csv'
    
    #Save new mouse and human data to file
    dfMouseTransformed.to_csv(os.path.join(outdir,fileMouseTx), index = False)
    dfHumanTransformed.to_csv(os.path.join(outdir,fileHumanTx), index = False)
    
    
    if args['voxeltransform'] == 'true':
        
        dfMouseVoxelTransformed = pd.DataFrame(net.predict_proba(X))
        dfMouseVoxelTransformed['Region'] = dfLabels[labelcol]
        
        file_voxel_human = ("HumanExpressionMatrix_"
                            "samples_pipeline_v1_labelled.csv")
        filepath_voxel_human = os.path.join(datadir, file_voxel_human)
        dfExprVoxelHuman = pd.read_csv(filepath_voxel_human)
        indLabelsHuman = dfExprVoxelHuman.columns.str.match('Region')
        dfInputVoxelHuman = dfExprVoxelHuman.loc[:, ~indLabels]
        X_VoxelHuman = dftx.fit_transform(dfInputVoxelHuman)['X']
        dfHumanVoxelTransformed = pd.DataFrame(net.predict_proba(X_VoxelHuman))
        dfHumanVoxelTransformed['Region'] = dfExprVoxelHuman[args['humandata'].capitalize()]
        
        fileMouseVoxelTx = "MLP_"+args['labels'].capitalize()+\
        "_Layers3"+\
        "_Units"+str(args['nunits'])+\
        "_L2"+str(args['L2'])+\
        "_MouseVoxelTx.csv"
        
        fileHumanVoxelTx = "MLP_"+args['labels'].capitalize()+\
        "_Layers3"+\
        "_Units"+str(args['nunits'])+\
        "_L2"+str(args['L2'])+\
        "_HumanVoxelTx.csv"

        dfMouseVoxelTransformed.to_csv(os.path.join(outdir,fileMouseVoxelTx), 
                                       index = False)
        dfHumanVoxelTransformed.to_csv(os.path.join(outdir,fileHumanVoxelTx), 
                                       index = False)
        
if __name__ == "__main__":
    main()

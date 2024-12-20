# Spatial-DC: a robust deep learning-based method for deconvolution of spatial proteomics

The framework of Spatial-DC, is described below.

<p align="center">
  <img width="70%" src=framework.jpg>
</p>

<p align="center"><strong>The framework of Spatial-DC. </strong></p>


## Stage I: Training the Deep Neural Network-Based Distribution Model (DIS)

In Stage I, Spatial-DC requires reference proteomics data **R** (e.g., cell-type proteomics) as input and obtains synthetic cell-type composition **Q** from a Dirichlet distribution (Methods). **R** and **Q** are then used to generate synthetic spots **M** and synthetic cell-type proteomic profile weights **N**. Subsequently, **Q** and **N** are stacked to produce synthetic labels **L**. **M** and **L** are used as training data and labels, respectively, to train a deep neural network-based distribution model **DIS**.

## Stage II: Training the Self-Supervised Model (SSM)

In Stage II, Spatial-DC utilizes spatial proteomics data, including spatial proteomic profiles **SEXP** and associated spatial coordinates **SCOO** as inputs. **SEXP** are fed into the trained **DIS** model to produce intermediate prediction **P**. Subsequently, **SEXP**, **SCOO**, and **P** are utilized to train a self-supervised model **SSM**, which is composed of:

- A variation graph autoencoder (VGAE) with encoder **G-E** and decoder **G-D**
- A deep autoencoder (DAE) with encoder **D-E** and decoder **D-D**
- A classifier (CLS)

Once trained, the model produces refined prediction **Y** from the classifier **CLS**, which are then split into cell-type composition **D** and cell-type proteomic profile weights **W**. By using **W**, along with the spatial proteomic profiles **SEXP**, cell-type proteomic profiles **E** of each spot are obtained.

## Abbreviations

- **c**: The number of cell types from reference data **R**
- **p**: The number of proteins shared between reference data **R** and spatial proteomic profiles **SEXP**
- **m**: The number of synthetic spots
- **n**: The number of spots from spatial proteomics **SEXP**
- **D-E/D-D**: Represents the encoder and decoder architecture of the autoencoder
- **G-E/G-D**: Represents the encoder and decoder architecture of the VGAE
```

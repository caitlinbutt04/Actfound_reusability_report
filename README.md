# ActFound reusability report
This repository contains the code for Reusability Report: evaluating the performance of a meta-learning foundation model on predicting the antibacterial activity of natural products. 
It is a fork of Actfound_demo (see below for original Actfound_demo README).

## Colabs for reusability report
In this reusability report, we fine-tuned five models with an antibacterial natural products dataset. The following colabs contain the code to fine-tune the models.

- Colab for Actfound
  
    - https://colab.research.google.com/drive/1c-yeku2Z9ZDQqtdyD2D0L3T4j_FBAova?usp=sharing
  
- Colab for Actfound with KNN-MAML and fusion method
  
    - https://colab.research.google.com/drive/1rsyZCtmOQ3OV2G2IIlMHysmCwhjh8L9M?usp=sharing
      
- Colab for MAML
  
    - https://colab.research.google.com/drive/17CLTovWtj7gJrMXzHhP3GqGUrXqrQ58m?usp=sharing
      
- Colab for ProtoNet
  
    - https://colab.research.google.com/drive/1KskTY91c_H4LE9mkSjkGMOC5oEknkef0?usp=sharing
      
- Colab for TransferQSAR
  
    - https://colab.research.google.com/drive/1gJcmdKbzzfMI4EykbpWhqNAeT2mpCylM?usp=sharing
  
# Actfound_demo

## Colab programmatic tool
⭐ We have now added a ready-to-use online programmatic tool on Colab.

😊 You can easily use Colab tools to fine-tune Actfound with few measured compounds and then use it to give prediction on your unmeasured compounds.

❤️ We also provide a metric in this colab to help you predict if Actfound can works well on your data with only few measured compounds for fine-tuning. (See Figure 3.g.h of our paper for more detailes)

- Colab for Actfound

    - https://colab.research.google.com/drive/1zn3U3xwLXZjZQYGVgtxCiEwt1hkho62l?usp=sharing 

- Colab for Actfound with KNN-MAML and fusion method

    - https://colab.research.google.com/drive/1eLWidAOWUqSCEcm0qM0Pf4IO1Ex7Ceal?usp=sharing

- Colab for finetuning full parameters of Actfound (can be more suitable when you have many measured compounds for finetuning, e.g. more than 300)

    - https://colab.research.google.com/drive/1EPx2tMHIdvhPbY8GyvuICVwVgwd3ikfJ?usp=sharing

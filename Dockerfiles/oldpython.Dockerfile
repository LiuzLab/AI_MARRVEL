FROM python:3.8

RUN pip install scikit-learn==1.1.2
RUN pip install pandas
RUN pip install scipy
RUN pip install xgboost


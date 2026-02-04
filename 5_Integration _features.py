import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# 读取数据
df = pd.read_csv('kernel-allfeatrue-filter.txt', sep='\t', index_col=0)
#df = pd.read_csv('merged_results.txt', sep='\t', index_col=0)

# 预处理函数
def preprocess_data(df):
    # 处理缺失值：用中位数填充
    df = df.fillna(df.median())
    
    # 处理P值类特征：转换为 -log10(p)，避免除零
    p_columns = [col for col in df.columns if ('-P' in col or 'eqtl' in col or '_pvalue' in col)]
    for col in p_columns:
        min_p = df[col][df[col] > 0].min()  # 找最小非零P值
        df[col] = df[col].replace(0, min_p/2)  # 将0替换为最小值的半值
        df[col] = -np.log10(df[col])
    
    # 处理TajimaD（取绝对值）
    if 'TajimaD' in df.columns:
        df['TajimaD'] = np.abs(df['TajimaD'])
    
    # 处理共表达系数（取绝对值）
    r_columns = [col for col in df.columns if 'Cor' in col]
    for col in r_columns:
        df[col] = np.abs(df[col])
    
    return df

# 标准化处理
def standardize_data(df):
    features = df.columns[1:]  # 假设第一列为基因名
    scaler = StandardScaler()
    df[features] = scaler.fit_transform(df[features])
    return df

# 综合评分方法
def calculate_scores(df):
    # 方法1：等权重求和
    features = df.columns[1:]
    df['Composite_Score'] = df[features].mean(axis=1)
    
    # 方法2：PCA主成分得分（推荐）
    pca = PCA(n_components=1)
    df['PCA_Score'] = pca.fit_transform(df[features])
    
    return df.sort_values(by='PCA_Score', ascending=False)

# 主流程
processed_df = preprocess_data(df.copy())
standardized_df = standardize_data(processed_df)
result_df = calculate_scores(standardized_df)

# 保存结果
result_df.to_csv('ranked_genes-allfeatrue-filter.txt', index=True, sep='\t')

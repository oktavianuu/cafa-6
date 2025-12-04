CAFA-6

# Project Name

## Dataset Setup

The full dataset is hosted on Kaggle. To use this repository:

### Option 1: Download from Kaggle directly
1. Visit: [Kaggle Dataset Link](https://www.kaggle.com/competitions/cafa-6-protein-function-prediction/data)
2. Download the dataset
3. Extract to the `data/` folder

### Option 2: Use Kaggle API
```bash
# Install kaggle API
pip install kaggle

# Configure with your API token
# (Get token from Kaggle Account Settings → API)
mkdir -p ~/.kaggle
cp kaggle.json ~/.kaggle/
chmod 600 ~/.kaggle/kaggle.json

# Download dataset
kaggle datasets download -d username/dataset-name -p data/
unzip data/dataset-name.zip -d data/

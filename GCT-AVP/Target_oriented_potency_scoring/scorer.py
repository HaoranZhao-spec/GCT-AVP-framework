import os
import pandas as pd
import argparse
import pathlib
import torch
from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained, MSATransformer
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset, random_split
import pickle
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from tqdm import tqdm

def to_fasta(from_csv_path, to_fasta_path):
    df = pd.read_csv(from_csv_path)
    seq_list = df['Sequence'].tolist()
    number = len(seq_list) + 1
    name = [i for i in range(1, number)]

    with open(to_fasta_path, 'w') as fasta_file:
        for i in name:
            fasta_lines = '>' + str(i)
            seq_lines = seq_list[i - 1]
            fasta_file.write(fasta_lines + '\n' + seq_lines + '\n')

def get_embedding(model_location, fasta_file, output_dir, toks_per_batch=4096, truncation_seq_length=1022, include='mean', repr_layers=[-1]):
    assert os.path.exists(fasta_file)
    fasta_file = pathlib.Path(fasta_file)
    output_dir = pathlib.Path(output_dir)
    model, alphabet = pretrained.load_model_and_alphabet(model_location)
    model.eval()
    if isinstance(model, MSATransformer):
        raise ValueError("This script currently does not handle models with MSA input (MSA Transformer).")
    if torch.cuda.is_available():
        model = model.cuda()
        print("Transferred model to GPU")

    dataset = FastaBatchedDataset.from_file(fasta_file)
    batches = dataset.get_batch_indices(toks_per_batch, extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(dataset, collate_fn=alphabet.get_batch_converter(truncation_seq_length), batch_sampler=batches)
    print(f"Read {fasta_file} with {len(dataset)} sequences")

    output_dir.mkdir(parents=True, exist_ok=True)
    return_contacts = "contacts" in include

    assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in repr_layers)
    repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in repr_layers]

    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(data_loader):
            print(f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)")
            if torch.cuda.is_available():
                toks = toks.to(device="cuda", non_blocking=True)

            out = model(toks, repr_layers=repr_layers, return_contacts=return_contacts)

            logits = out["logits"].to(device="cpu")
            representations = {layer: t.to(device="cpu") for layer, t in out["representations"].items()}
            if return_contacts:
                contacts = out["contacts"].to(device="cpu")

            for i, label in enumerate(labels):
                output_file = output_dir / f"{label}.pt"
                output_file.parent.mkdir(parents=True, exist_ok=True)
                result = {"label": label}
                truncate_len = min(truncation_seq_length, len(strs[i]))
                if "per_tok" in include:
                    result["representations"] = {layer: t[i, 1 : truncate_len + 1].clone() for layer, t in representations.items()}
                if "mean" in include:
                    result["mean_representations"] = {layer: t[i, 1 : truncate_len + 1].mean(0).clone() for layer, t in representations.items()}
                if "bos" in include:
                    result["bos_representations"] = {layer: t[i, 0].clone() for layer, t in representations.items()}
                if return_contacts:
                    result["contacts"] = contacts[i, : truncate_len, : truncate_len].clone()

                torch.save(result, output_file)

def load_embeding(from_folder_path, from_csv_path, to_device):
    folder_path = from_folder_path
    csv_file_path = from_csv_path
    device = torch.device('cuda' if torch.cuda.is_available() and to_device == 'cuda' else 'cpu')
    print(f"Running on {device}")

    file_names = sorted([f for f in os.listdir(folder_path) if f.endswith('.pt')], key=lambda x: int(x.split('.')[0]))
    mean_representations_list = []

    for file_name in file_names:
        file_path = os.path.join(folder_path, file_name)
        mean_representation = torch.load(file_path,  map_location=device, weights_only=True)['mean_representations'][36].cpu().numpy().tolist()
        mean_representations_list.append(mean_representation)

    mean_representations_df = pd.DataFrame(mean_representations_list)
    additional_data_df = pd.read_csv(csv_file_path)

    if len(mean_representations_df) != len(additional_data_df):
        raise ValueError("The number of rows in the additional data CSV does not match the number of .pt files")

    final_df = pd.concat([additional_data_df['Sequence'], mean_representations_df], axis=1)

    if 'EC50' in additional_data_df.columns:
        final_df['EC50] = additional_data_df['EC50']
    return final_df

class LSTMModel(nn.Module):
    def __init__(self, input_size, hidden_size, num_layers, output_size, dropout_rate):
        super(LSTMModel, self).__init__()
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers, batch_first=True, dropout=dropout_rate)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc = nn.Linear(hidden_size, output_size)
    
    def forward(self, x):
        h_0 = torch.zeros(self.lstm.num_layers, x.size(0), self.lstm.hidden_size).to(x.device)
        c_0 = torch.zeros(self.lstm.num_layers, x.size(0), self.lstm.hidden_size).to(x.device)
        out, _ = self.lstm(x, (h_0, c_0))
        out = self.dropout(out[:, -1, :])
        out = self.fc(out)
        return out

def train_lstm_model(
    train_data_df, 
    val_data_df, 
    scaler_save_path, 
    model_save_path, 
    to_device,
    hidden_size=128,
    num_layers=2,
    dropout_rate=0.7,
    batch_size=64,
    epochs=100,
    lr=1e-3,
    patience=10 
):
    """
    训练LSTM模型并保存标度器、模型
    :param train_data_df: 训练集（含Sequence、嵌入列、EC50列）
    :param val_data_df: 验证集（结构同训练集）
    :param scaler_save_path: 标度器保存路径（如./scaler.pkl）
    :param model_save_path: 模型保存路径（如./lstm_model.pt）
    :param to_device: 训练设备（cuda/cpu）
    """
    device = torch.device('cuda' if torch.cuda.is_available() and to_device == 'cuda' else 'cpu')
    print(f"Training on {device}")

    X_train = train_data_df.iloc[:, 1:-1].values  # 嵌入列（排除Sequence和EC50）
    y_train = train_data_df['EC50'].values
    X_val = val_data_df.iloc[:, 1:-1].values
    y_val = val_data_df['EC50'].values

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_val_scaled = scaler.transform(X_val)
    with open(scaler_save_path, 'wb') as f:
        pickle.dump(scaler, f)
    print(f"Scaler saved to {scaler_save_path}")

    X_train_tensor = torch.tensor(X_train_scaled, dtype=torch.float32).unsqueeze(1)  # [batch, 1, input_size]
    y_train_tensor = torch.tensor(y_train, dtype=torch.float32).unsqueeze(1)
    X_val_tensor = torch.tensor(X_val_scaled, dtype=torch.float32).unsqueeze(1)
    y_val_tensor = torch.tensor(y_val, dtype=torch.float32).unsqueeze(1)

    train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
    val_dataset = TensorDataset(X_val_tensor, y_val_tensor)

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

    input_size = X_train_scaled.shape[1]
    model = LSTMModel(
        input_size=input_size,
        hidden_size=hidden_size,
        num_layers=num_layers,
        output_size=1,
        dropout_rate=dropout_rate
    ).to(device)

    criterion = nn.MSELoss()  
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=5, factor=0.5)

    best_val_loss = float('inf')
    early_stop_count = 0
    train_losses = []
    val_losses = []

    for epoch in range(epochs):

        model.train()
        train_loss = 0.0
        for X_batch, y_batch in tqdm(train_loader, desc=f"Epoch {epoch+1}/{epochs}"):
            X_batch = X_batch.to(device)
            y_batch = y_batch.to(device)
            
            optimizer.zero_grad()
            outputs = model(X_batch)
            loss = criterion(outputs, y_batch)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item() * X_batch.size(0)
        
        train_loss /= len(train_loader.dataset)
        train_losses.append(train_loss)

        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch = X_batch.to(device)
                y_batch = y_batch.to(device)
                outputs = model(X_batch)
                loss = criterion(outputs, y_batch)
                val_loss += loss.item() * X_batch.size(0)
        
        val_loss /= len(val_loader.dataset)
        val_losses.append(val_loss)

        scheduler.step(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            early_stop_count = 0

            torch.save(model.state_dict(), model_save_path)
            print(f"Epoch {epoch+1}: Val loss improved to {best_val_loss:.4f}, saving model")
        else:
            early_stop_count += 1
            if early_stop_count >= patience:
                print(f"Early stopping at epoch {epoch+1} (no improvement for {patience} epochs)")
                break

        print(f"Epoch {epoch+1} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

    plt.figure(figsize=(10, 6))
    plt.plot(train_losses, label='Train Loss')
    plt.plot(val_losses, label='Val Loss')
    plt.xlabel('Epoch')
    plt.ylabel('MSE Loss')
    plt.title('Training & Validation Loss')
    plt.legend()
    plt.savefig('./train_val_loss.png')
    plt.show()

    print(f"Best model saved to {model_save_path}")
    return model

def get_predicted_EC50(new_data, from_csv_path, scaler_data_path, model_path, result_path, to_device):
    device = torch.device('cuda' if torch.cuda.is_available() and to_device == 'cuda' else 'cpu')
    print(f"Running on {device}")

    X_new = new_data.iloc[:, 1:].values
    with open(scaler_data_path, 'rb') as f:
        scaler = pickle.load(f)
    X_new = scaler.transform(X_new)

    X_new = torch.tensor(X_new, dtype=torch.float32)

    new_dataset = TensorDataset(X_new.unsqueeze(1))
    new_loader = DataLoader(new_dataset, batch_size=64, shuffle=False)

    model = LSTMModel(input_size=X_new.size(1), hidden_size=128, num_layers=2, output_size=1, dropout_rate=0.7).to(device)
    model_state_dict = torch.load(model_path, map_location=device, weights_only=True)
    model.load_state_dict(model_state_dict)
    model.eval()

    new_predictions = []
    with torch.no_grad():
        for X_batch in new_loader:
            X_batch = X_batch[0].to(device)
            outputs = model(X_batch)
            new_predictions.extend(outputs.squeeze().cpu().numpy())

    new_predictions = pd.DataFrame(new_predictions, columns=['Predicted Values'])
    df_info  = pd.read_csv(from_csv_path)
    df_merged = pd.concat([df_info, new_predictions], axis=1)

    df_merged.to_csv(result_path, index=False)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--mode', required=True, type=str, choices=['train', 'predict'], help='Mode: train/predict')
    parser.add_argument('--from_csv_path', required=True, type=str, help='Path to input CSV (train: with EC50; predict: only Sequence)')
    parser.add_argument('--to_fasta_path', type=str, default='./data/seqs.fasta')
    parser.add_argument('--esm_model_location', type=str, default='esm2_t36_3B_UR50D')
    parser.add_argument('--output_dir', type=str, default='./data/esm_output/')
    parser.add_argument('--repr_layers', type=int, default=36)
    parser.add_argument('--to_device', default='cuda', type=str)

    parser.add_argument('--scaler_save_path', type=str, default='./data/scaler.pkl', help='Path to save scaler (train mode)')
    parser.add_argument('--model_save_path', type=str, default='./data/lstm_model.pt', help='Path to save LSTM model (train mode)')
    parser.add_argument('--val_split', type=float, default=0.2, help='Validation split ratio (train mode)')
    parser.add_argument('--epochs', type=int, default=100, help='Training epochs (train mode)')
    parser.add_argument('--batch_size', type=int, default=64, help='Batch size (train mode)')
    parser.add_argument('--lr', type=float, default=1e-3, help='Learning rate (train mode)')

    parser.add_argument('--scaler_data_path', type=str, help='Path to scaler (predict mode)')
    parser.add_argument('--model_path', type=str, help='Path to trained LSTM model (predict mode)')
    parser.add_argument('--result_path', type=str, help='Path to save prediction result (predict mode)')

    args = parser.parse_args()

    print(f'[i] to_fasta')
    to_fasta(args.from_csv_path, args.to_fasta_path)
    print(f'[i] get_embedding')
    get_embedding(args.esm_model_location, args.to_fasta_path, args.output_dir, repr_layers=[args.repr_layers])
    print(f'[i] load_embeding')
    df = load_embeding(args.output_dir, args.from_csv_path, args.to_device)

    if args.mode == 'train':

        if 'EC50' not in df.columns:
            raise ValueError("Training CSV must contain 'EC50' column as label")
        
        train_df, val_df = train_test_split(df, test_size=args.val_split, random_state=42)
        print(f"Train set: {len(train_df)} samples | Val set: {len(val_df)} samples")
        
        print(f'[i] train_lstm_model')
        train_lstm_model(
            train_data_df=train_df,
            val_data_df=val_df,
            scaler_save_path=args.scaler_save_path,
            model_save_path=args.model_save_path,
            to_device=args.to_device,
            epochs=args.epochs,
            batch_size=args.batch_size,
            lr=args.lr
        )

    elif args.mode == 'predict':

        if not all([args.scaler_data_path, args.model_path, args.result_path]):
            raise ValueError("Predict mode requires --scaler_data_path, --model_path, --result_path")
        
        print(f'[i] get_predicted_EC50')
        get_predicted_EC50(df, args.from_csv_path, args.scaler_data_path, args.model_path, args.result_path, args.to_device)

if __name__ == '__main__':
    main()
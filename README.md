# Tugas Paper GSTAR - Mata Kuliah Teknik Peramalan (T)

## Anggota Kelompok

| Nama | NRP |
| --- | --- | 
| Abdullah Sultan Barizy | 5025241092 | 
| Farikh Muhammad Fauzan | 5025241135 |
| Willy Dava Nugraha | 5025241090 |

## Deskripsi Dataset 
- Dataset diambil dari `https://dataonline.bmkg.go.id/` dengan rentang data 12 bulan (Mei 2025 - April 2026) dari tiga stasiun cuaca dari dua kota, yaitu dari Kota Surabaya (Stasiun Meteorologi Perak I dan Stasiun Meteorologi Maritim Tanjung Perak) dan Kabupaten Sidoarjo (Stasiun Meteorologi Juanda).
- Dataset mentah dapat dilihat pada folder `data/raw`. Setiap file data mentah memiliki keterbatasan pengambilan hanya berisi rentang maksimal 30 hari dari satu stasiun cuaca.
- Dataset yang telah diproses pada rentang 12 bulan ((Mei 2025 - April 2026) dan 2 bulan (Maret 2026 - April 2026) dapat dilihat pada folder `data/processed`.
- Plot dataset yang sudah diproses selama 12 bulan ditampilkan pada grafik berikut:
  <img width="1320" height="585" alt="image" src="https://github.com/user-attachments/assets/8edcf85f-8430-4c23-8cc2-be6276e90b7c" />

## Uji Stasioner 
- Uji stasioner dengan menggunakan CADF dengan kode pada path `scripts/01_Stationarity_CADF_optimal.R` dengan pengujian sampai lag ke 10 diperoleh stasioner sampai lag ke-4.
  <img width="606" height="126" alt="image" src="https://github.com/user-attachments/assets/122a576c-3397-4e1f-b780-6a87afc611bc" />

## Plot ACF dan PACF
<img width="655" height="545" alt="image" src="https://github.com/user-attachments/assets/a3585828-c5ba-4c21-9bb2-8c491e7bf8e2" />

## Summary Modeling
<img width="441" height="318" alt="image" src="https://github.com/user-attachments/assets/cf32fd12-da5b-4ab7-aa6a-9bfeeda1098b" />

## Performance Training
- Evaluasi dengan MSE dan MAPE
<img width="342" height="194" alt="image" src="https://github.com/user-attachments/assets/5fa472f8-7084-4bb1-b5e9-f978f3bd2f24" />
  
- Evaluasi dengan $$R^2$$
<img width="606" height="336" alt="image" src="https://github.com/user-attachments/assets/6cc33443-f5a7-4924-be0a-c94781bdc378" />

## Performance Testing
- Evaluasi dengan MSE dan MAPE
<img width="346" height="388" alt="image" src="https://github.com/user-attachments/assets/34722a71-943d-491c-89f4-fb4492bfde3a" />

- Evaluasi dengan $$R^2$$
<img width="477" height="259" alt="image" src="https://github.com/user-attachments/assets/f0bbfe34-1e91-4053-a3cf-ea68d6e7f9fa" />

## Plot Model
<img width="1331" height="545" alt="image" src="https://github.com/user-attachments/assets/c149f791-ca16-4b20-bf6d-898156b32b36" />





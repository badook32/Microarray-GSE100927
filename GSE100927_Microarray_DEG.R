# 필요한 라이브러리 로드
library(oligo)        # oligo 패키지 사용
library(limma)
library(hgu133plus2.db)
library(pheatmap)
library(ggplot2)

# Step 1: 데이터 로드
setwd("~/Desktop/setting/기존 데이터셋/GSE100927_RAW(병변:동맥 비교)")  # 데이터 디렉토리 설정

# 1. 파일 목록 읽어오기
files <- list.files(path = "~/Desktop/setting/기존 데이터셋/GSE100927_RAW(병변:동맥 비교)", pattern = "*.txt.gz", full.names = TRUE)

# 파일 목록 확인
print(files)

# 샘플 정보 매칭 (샘플 정보 데이터)
sample_info <- data.frame(
  GSM = c("GSM2696609", "GSM2696610", "GSM2696611", "GSM2696612", "GSM2696613", "GSM2696614", 
          "GSM2696615", "GSM2696616", "GSM2696617", "GSM2696618", "GSM2696619", "GSM2696620",
          "GSM2696621", "GSM2696622", "GSM2696623", "GSM2696624", "GSM2696625", "GSM2696626", 
          "GSM2696627", "GSM2696628", "GSM2696629", "GSM2696630", "GSM2696631", "GSM2696632",
          "GSM2696633", "GSM2696634", "GSM2696635", "GSM2696636", "GSM2696637", "GSM2696638",
          "GSM2696639", "GSM2696640", "GSM2696641", "GSM2696642", "GSM2696643", "GSM2696644",
          "GSM2696645", "GSM2696646", "GSM2696647", "GSM2696648", "GSM2696649", "GSM2696650",
          "GSM2696651", "GSM2696652", "GSM2696653", "GSM2696654", "GSM2696655", "GSM2696656",
          "GSM2696657", "GSM2696658", "GSM2696659", "GSM2696660", "GSM2696661", "GSM2696662",
          "GSM2696663", "GSM2696664", "GSM2696665", "GSM2696666", "GSM2696667", "GSM2696668",
          "GSM2696669", "GSM2696670", "GSM2696671", "GSM2696672", "GSM2696673", "GSM2696674",
          "GSM2696675", "GSM2696676", "GSM2696677", "GSM2696678", "GSM2696679", "GSM2696680",
          "GSM2696681", "GSM2696682", "GSM2696683", "GSM2696684", "GSM2696685", "GSM2696686",
          "GSM2696687", "GSM2696688", "GSM2696689", "GSM2696690", "GSM2696691", "GSM2696692",
          "GSM2696693", "GSM2696694", "GSM2696695", "GSM2696696", "GSM2696697", "GSM2696698",
          "GSM2696699", "GSM2696700", "GSM2696701", "GSM2696702", "GSM2696703", "GSM2696704",
          "GSM2696705", "GSM2696706", "GSM2696707", "GSM2696708", "GSM2696709", "GSM2696710",
          "GSM2696711", "GSM2696712"),
  Sample_Name = c("H1_87: Femoral artery_3", "H10_109: Femoral artery_7", "H102_8: Femoral artery_Control 26", 
                  "H103_90: Carotid artery_Control 14", "H104_154: Infra-popliteal artery_Control 01", 
                  "H105_25: Carotid artery_91", "H106_89b: Carotid artery_Control 10", 
                  "H108_34: Infra-popliteal artery_Control 16", "H109_78: Carotid artery_Control 19", 
                  "H110_54: Carotid artery_Control 12", "H111_40: Infra-popliteal artery_1", 
                  "H112_26: Carotid artery_10", "H113_77: Carotid artery_Control 27", 
                  "H114_84: Carotid artery_84", "H116_49: Carotid artery_2", "H117_56: Femoral artery_Control 08", 
                  "H121_92: Femoral artery_Control 21", "H122_169: Femoral artery_11", 
                  "H123_1: Carotid artery_Control 22", "H126_146: Carotid artery_27", 
                  "H128_73: Femoral artery_25", "H129_168: Carotid artery_82", 
                  "H13_193: Infra-popliteal artery_10", "H130_150: Infra-popliteal artery_19", 
                  "H131_28: Femoral artery_1", "H132_88: Femoral artery_31", "H133_44b: Femoral artery_Control 05", 
                  "H135_17: Femoral artery_Control 01", "H136_99: Carotid artery_55", 
                  "H138_185b: Femoral artery_27", "H139_13: Infra-popliteal artery_12", 
                  "H14_23: Femoral artery_9", "H140_190: Femoral artery_23", 
                  "H141_93: Infra-popliteal artery_Control 19", "H142_189: Carotid artery_32", 
                  "H143_39: Femoral artery_18", "H144_81: Carotid artery_44", 
                  "H145_46: Carotid artery_40", "H146_15: Carotid artery_Control 05", 
                  "H147_112: Infra-popliteal artery_Control 25", "H148_104: Carotid artery_95", 
                  "H149_174: Infra-popliteal artery_Control 18", "H154_194: Infra-popliteal artery_Control 11", 
                  "H155_103: Carotid artery_33", "H156_36: Carotid artery_12", 
                  "H157_184: Femoral artery_40", "H159_140: Femoral artery_28", 
                  "H160_172: Femoral artery_41G", "H162_59: Carotid artery_70", 
                  "H167_21: Carotid artery_42", "H168_182: Femoral artery_4", 
                  "H17_200: Femoral artery_Control 25", "H173_170: Femoral artery_16", 
                  "H177_153: Femoral artery_Control 22", "H179_158: Carotid artery_22", 
                  "H18_75: Infra-popliteal artery_20", "H182_185: Femoral artery_27", 
                  "H184_196: Infra-popliteal artery_17", "H185_163: Carotid artery_6", 
                  "H188_149: Femoral artery_20", "H190_44: Femoral artery_Control 05", 
                  "H191_82: Carotid artery_29", "H194_201: Infra-popliteal artery_Control 11", 
                  "H21_53: Infra-popliteal artery_3", "H22_98: Carotid artery_92", 
                  "H23_9: Carotid artery_5", "H26_47: Carotid artery_26", 
                  "H27_183: Femoral artery_32", "H28_61: Femoral artery_17", 
                  "H29_33: Femoral artery_Control 11", "H32_10: Carotid artery_81", 
                  "H34_26b: Carotid artery_10", "H35_97: Carotid artery_38", 
                  "H36_173: Infra-popliteal artery_Control 08", "H38_76: Infra-popliteal artery_9", 
                  "H40_24: Infra-popliteal artery_17", "H51_151: Infra-popliteal artery_7", 
                  "H52_141: Carotid artery_44", "H53_52: Femoral artery_15", 
                  "H55_63: Carotid artery_Control 21", "H58_102: Carotid artery_21", 
                  "H59_192: Infra-popliteal artery_5", "H6_4: Femoral artery_Control 25", 
                  "H60_186: Femoral artery_2", "H61_155: Infra-popliteal artery_Control 10", 
                  "H65_94: Infra-popliteal artery_Control 23", "H66_91: Femoral artery_Control 06", 
                  "H69_20: Carotid artery_37", "H72_198: Carotid artery_Control 22", 
                  "H73_29: Infra-popliteal artery_11", "H74_65: Femoral artery_Control 17", 
                  "H79_89: Carotid artery_Control 10", "H8_69: Carotid artery_16", 
                  "H84_12: Femoral artery_14", "H87_171: Femoral artery_22", 
                  "H88_110: Femoral artery_43", "H89_74: Femoral artery_29", 
                  "H90_62: Infra-popliteal artery_4", "H91_111: Femoral artery_12", 
                  "H93_96: Carotid artery_48", "H96_31: Carotid artery_Control 16", 
                  "H97_18: Infra-popliteal artery_Control 12", "H98_42: Carotid artery_Control 01", 
                  "H99_191: Infra-popliteal artery_8")
)

# 2. 파일 이름 매칭 및 수정
file_names <- list.files(path = "~/Desktop/setting/기존 데이터셋/GSE100927_RAW(병변:동맥 비교)", 
                         pattern = "*.txt.gz", full.names = TRUE)

# 파일 이름 수정 (GSM 번호 매칭)
new_file_names <- sapply(file_names, function(x) {
  gsm_code <- gsub(".*/(GSM\\d+)_.*", "\\1", x)  # GSM 번호 추출
  matched_sample <- sample_info$Sample_Name[sample_info$GSM == gsm_code]  # 매칭되는 샘플 이름
  new_name <- gsub(gsm_code, matched_sample, x)  # 기존 파일명에서 GSM 번호를 샘플 이름으로 변경
  return(new_name)
})

# 결과 확인
head(new_file_names)

# 파일 읽기 및  데이터 병합
library(data.table)

# "FEATURES" 섹션 이후 데이터를 읽고, 샘플 이름을 열 이름으로 지정
read_microarray_data <- function(file) {
  con <- gzfile(file, "r")
  lines <- readLines(con)
  close(con)
  
  feature_line <- grep("^FEATURES", lines)[1]
  if (is.na(feature_line)) return(NULL)
  
  # 데이터를 fread로 읽고, 필요한 열 선택
  data <- fread(text = readLines(gzfile(file)), skip = feature_line - 1, sep = "\t", fill = TRUE, data.table = TRUE)
  
  # 샘플 이름을 파일명에서 추출하여 열 이름으로 지정
  sample_name <- gsub(".*/|\\.txt\\.gz", "", file)  # 파일명에서 샘플 ID 추출
  colnames(data)[which(colnames(data) == "gProcessedSignal")] <- sample_name
  
  # ProbeName 기준으로 평균값만 남김
  data <- data[, .(MeanExpression = mean(get(sample_name), na.rm = TRUE)), by = ProbeName]
  setnames(data, "MeanExpression", sample_name)  # 열 이름 변경
  
  return(data)
}

# 모든 파일을 읽어 리스트로 저장
raw_data_list <- lapply(files, read_microarray_data)

# NA 값이 포함된 경우 제거
raw_data_list <- raw_data_list[!sapply(raw_data_list, is.null)]

# 중복 없이 병합 (ProbeName을 기준으로)
merged_data <- Reduce(function(x, y) merge(x, y, by = "ProbeName", all = TRUE), raw_data_list)

# 올바른 형태인지 확인
print(dim(merged_data))  # 행(Probe 개수) x 열(샘플 개수 + 1(ProbeName))
head(merged_data)


# 데이터를 데이터프레임으로 변환
normalized_data <- as.data.frame(merged_data)

# ProbeName을 rownames로 설정
rownames(normalized_data) <- normalized_data$ProbeName
normalized_data$ProbeName <- NULL  # 기존 열 제거

# Log2 변환
normalized_data <- log2(normalized_data + 1)

# 확인
print(dim(normalized_data))  # 샘플 개수가 여러 개 나와야 정상
head(normalized_data)


# pheno_data 생성
pheno_data <- data.frame(
  Sample = colnames(normalized_data),
  Condition = rep(c("Femoral artery", "Femoral artery_Control", 
                    "Carotid artery", "Carotid artery_Control",
                    "Infra-popliteal artery", "Infra-popliteal artery_Control"),
                  length.out = length(colnames(normalized_data))),
  row.names = colnames(normalized_data)
)

# 확인
print(dim(pheno_data))  # 정상적으로 샘플 수가 반영되었는지 확인
head(pheno_data)

print(colnames(normalized_data)[1:10])  # 앞 10개 샘플 이름 확인

# 1) normalized_data의 열 이름에서 GSM 코드 부분만 추출
extracted_gsm <- gsub("^(GSM\\d+).*", "\\1", colnames(normalized_data))

# 2) 추출한 GSM 코드가 sample_info에 존재하는지 확인
gsm_match_check <- extracted_gsm %in% sample_info$GSM
sum(gsm_match_check)  # ✅ TRUE 개수가 104개여야 정상

# 3)  NA가 포함된 경우 문제 있는 GSM 코드 출력
if (sum(!gsm_match_check) > 0) {
  print(extracted_gsm[!gsm_match_check])  # 매칭되지 않은 GSM 코드 출력
}

# GSM 코드와 샘플 이름을 매핑하는 named vector 생성
gsm_to_sample <- setNames(sample_info$Sample_Name, sample_info$GSM)

# 정상적으로 생성되었는지 확인
print(gsm_to_sample[1:5])  # 앞 5개 샘플만 확인

# 매칭이 잘 되는 경우, 올바른 샘플 이름 적용
colnames(normalized_data) <- gsm_to_sample[extracted_gsm]

# 변경 후 NA 개수 확인 (0이어야 정상)
sum(is.na(colnames(normalized_data)))

pheno_data <- data.frame(
  Sample = colnames(normalized_data),
  Condition = rep(c("Femoral artery", "Femoral artery_Control", 
                    "Carotid artery", "Carotid artery_Control",
                    "Infra-popliteal artery", "Infra-popliteal artery_Control"),
                  length.out = length(colnames(normalized_data))),
  row.names = colnames(normalized_data)
)

# 확인
print(dim(pheno_data))  # 샘플 수가 정상인지 확인
head(pheno_data)



# 1) 기관 이름만 남기기 (번호 및 불필요한 부분 제거)
pheno_data$Condition <- gsub(".*: (.+)", "\\1", pheno_data$Sample)  # 앞부분 제거
pheno_data$Condition <- gsub("_[0-9]+.*$", "", pheno_data$Condition)  # 숫자 뒤 전체 삭제
pheno_data$Condition <- gsub("\\s+$", "", pheno_data$Condition)  # 공백 제거

# 2) Control 샘플 변환 (뒤쪽 숫자/문자 삭제하여 "_Control"로 통일)
pheno_data$Condition <- ifelse(grepl("Control", pheno_data$Sample),
                               gsub("\\s*Control.*$", "_Control", pheno_data$Condition),  
                               pheno_data$Condition)

# 변환 후 최종 확인
table(pheno_data$Condition)

# __Control 되면 이 코드 실행
# 1) 기관 이름만 남기기 (번호 및 불필요한 부분 제거)
pheno_data$Condition <- gsub(".*: (.+)", "\\1", pheno_data$Sample)  # 앞부분 제거
pheno_data$Condition <- gsub("_[0-9]+.*$", "", pheno_data$Condition)  # 숫자 뒤 전체 삭제
pheno_data$Condition <- gsub("\\s+$", "", pheno_data$Condition)  # 공백 제거

# 2) Control 샘플 변환 (뒤쪽 숫자/문자 삭제하여 "_Control"로 통일)
pheno_data$Condition <- ifelse(grepl("Control", pheno_data$Sample),
                               gsub("\\s*Control.*$", "", pheno_data$Condition),  # "Control" 이후 부분 삭제
                               pheno_data$Condition)
pheno_data$Condition <- ifelse(grepl("Control", pheno_data$Sample),
                               paste0(pheno_data$Condition, "_Control"),  # "_Control" 추가
                               pheno_data$Condition)

# 3) 잘못된 "__Control" 수정 → "_Control"로 통일
pheno_data$Condition <- gsub("__Control", "_Control", pheno_data$Condition)

# 변환 후 최종 확인
table(pheno_data$Condition)





# A. 대퇴동맥 (Femoral artery) vs 대조군 (Control)
# 대퇴동맥(Femoral artery) vs 대조군(Femoral artery_Control) 비교
femoral_vs_con <- subset(pheno_data, grepl("Femoral artery", Sample))
femoral_vs_con$Condition <- factor(femoral_vs_con$Condition, levels = c("Femoral artery", "Femoral artery_Control"))

# Limma 분석 수행
design_femoral <- model.matrix(~0 + Condition, data = femoral_vs_con)
colnames(design_femoral) <- c("Femoral", "Control")
fit_femoral <- lmFit(normalized_data[, rownames(femoral_vs_con)], design_femoral)
contrast_femoral <- makeContrasts(Femoral - Control, levels = design_femoral)
fit_femoral2 <- contrasts.fit(fit_femoral, contrast_femoral)
fit_femoral2 <- eBayes(fit_femoral2)

# 결과 출력
topTable(fit_femoral2, adjust = "BH", p.value = 0.05)



# B. 경동맥 (Carotid artery) vs 대조군 (Control)
# 경동맥(Carotid artery) vs 대조군(Carotid artery_Control) 비교
carotid_vs_con <- subset(pheno_data, grepl("Carotid artery", Sample))
carotid_vs_con$Condition <- factor(carotid_vs_con$Condition, levels = c("Carotid artery", "Carotid artery_Control"))

# Limma 분석 수행
design_carotid <- model.matrix(~0 + Condition, data = carotid_vs_con)
colnames(design_carotid) <- c("Carotid", "Control")
fit_carotid <- lmFit(normalized_data[, rownames(carotid_vs_con)], design_carotid)
contrast_carotid <- makeContrasts(Carotid - Control, levels = design_carotid)
fit_carotid2 <- contrasts.fit(fit_carotid, contrast_carotid)
fit_carotid2 <- eBayes(fit_carotid2)

# 결과 출력
topTable(fit_carotid2, adjust = "BH", p.value = 0.05)



# 3. 하지동맥(Infra-popliteal artery) vs 대조군(Infra-popliteal artery_Control) 비교
infra_vs_con <- subset(pheno_data, grepl("Infra-popliteal artery", Sample))
infra_vs_con$Condition <- factor(infra_vs_con$Condition, levels = c("Infra-popliteal artery", "Infra-popliteal artery_Control"))

# Limma 분석 수행
design_infra <- model.matrix(~0 + Condition, data = infra_vs_con)

# R 변수명으로 변환 (자동 수정)
colnames(design_infra) <- make.names(colnames(design_infra))

# "Condition" 제거하여 올바른 변수명으로 수정
colnames(design_infra) <- gsub("^Condition", "", colnames(design_infra))

# 변수명 확인 (올바르게 변환되었는지)
print(colnames(design_infra))  # "Infra.popliteal.artery", "Infra.popliteal.artery_Control"

# 새로운 변수명을 반영하여 lmFit() 수행
fit_infra <- lmFit(normalized_data[, rownames(infra_vs_con)], design_infra)

# 올바르게 변환된 변수명을 사용하여 contrast 설정
contrast_infra <- makeContrasts(Infra.popliteal.artery - Infra.popliteal.artery_Control, levels = design_infra)
fit_infra2 <- contrasts.fit(fit_infra, contrast_infra)
fit_infra2 <- eBayes(fit_infra2)

# 최종 결과 출력
topTable(fit_infra2, adjust = "BH", p.value = 0.05)




# D. 모든 대조군 (All control) vs 나머지 샘플 (Other samples)
# Control 샘플 여부를 Condition에 반영
pheno_data$Condition <- ifelse(grepl("_Control$", pheno_data$Condition), "Control", "Other")

# 정상적으로 Control과 Other로 분류되었는지 확인
table(pheno_data$Condition)

# 모델 매트릭스 생성
design_all <- model.matrix(~0 + Condition, data = pheno_data)

# 자동 생성된 변수명 확인
print(colnames(design_all))  # "ConditionControl", "ConditionOther"여야 정상

# "Condition" 부분 제거하여 변수명 정리
colnames(design_all) <- gsub("^Condition", "", colnames(design_all))

# 새로운 변수명을 반영하여 lmFit() 수행
fit_all <- lmFit(normalized_data[, rownames(pheno_data)], design_all)

# 올바르게 변환된 변수명을 사용하여 contrast 설정
contrast_all <- makeContrasts(Other - Control, levels = design_all)
fit_all2 <- contrasts.fit(fit_all, contrast_all)
fit_all2 <- eBayes(fit_all2)

# 최종 결과 출력
topTable(fit_all2, adjust = "BH", p.value = 0.05)

# Limma 분석 수행 후 결과 저장
top_genes_femoral <- topTable(fit_femoral2, adjust = "BH", p.value = 0.05, number = Inf)
top_genes_carotid <- topTable(fit_carotid2, adjust = "BH", p.value = 0.05, number = Inf)
top_genes_infra <- topTable(fit_infra2, adjust = "BH", p.value = 0.05, number = Inf)
top_genes_all <- topTable(fit_all2, adjust = "BH", p.value = 0.05, number = Inf)


# 필요한 패키지
library(biomaRt)
library(dplyr)

# 1) Ensembl 마트 설정 (미러 서버 이용)
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      version = 108)

# 2) topTable 결과에서 probe ID 가져오기
top_genes_femoral$ID <- rownames(top_genes_femoral)
top_genes_carotid$ID <- rownames(top_genes_carotid)
top_genes_infra$ID <- rownames(top_genes_infra)
top_genes_all$ID <- rownames(top_genes_all)

# rownames 제거 (필요 시)
rownames(top_genes_femoral) <- NULL
rownames(top_genes_carotid) <- NULL
rownames(top_genes_infra) <- NULL
rownames(top_genes_all) <- NULL

# 3) 모든 probe ID 합치기 (중복 제거)
all_probe_ids <- unique(c(top_genes_femoral$ID,
                          top_genes_carotid$ID,
                          top_genes_infra$ID,
                          top_genes_all$ID))

# 4️) biomaRt로 매핑 수행
mapped <- getBM(
  attributes = c("agilent_sureprint_g3_ge_8x60k_v2", "hgnc_symbol", "entrezgene_id"),
  filters = "agilent_sureprint_g3_ge_8x60k_v2",
  values = all_probe_ids,
  mart = ensembl
)

# 5️) 결과 병합 (각 비교군에 매핑된 심볼 추가)
top_genes_femoral <- left_join(top_genes_femoral, mapped, by = c("ID" = "agilent_sureprint_g3_ge_8x60k_v2"))
top_genes_carotid <- left_join(top_genes_carotid, mapped, by = c("ID" = "agilent_sureprint_g3_ge_8x60k_v2"))
top_genes_infra   <- left_join(top_genes_infra, mapped, by = c("ID" = "agilent_sureprint_g3_ge_8x60k_v2"))
top_genes_all     <- left_join(top_genes_all, mapped, by = c("ID" = "agilent_sureprint_g3_ge_8x60k_v2"))

# 6) 확인
head(top_genes_femoral)
head(top_genes_carotid)
head(top_genes_infra)
head(top_genes_all)

# 7️) 필요한 경우 저장도 가능
# write.csv(top_genes_femoral, "top_genes_femoral_biomart.csv", row.names = FALSE)




library(ggplot2)
library(ggpubr)
library(dplyr)

# 1) OASL에 해당하는 probe ID 찾기
oasl_probe <- mapped$agilent_sureprint_g3_ge_8x60k_v2[mapped$hgnc_symbol == "OASL"]

# 2️) OASL 발현값 추출
oasl_expr <- as.numeric(normalized_data[oasl_probe, ])

# 3️) pheno_data에 OASL 발현 추가
pheno_data$OASL <- oasl_expr

# 4️) Condition을 CONTROL / AS로 단순화
pheno_data$Diagnosis <- ifelse(pheno_data$Condition == "Control", "CONTROL", "AS")
pheno_data$Diagnosis <- factor(pheno_data$Diagnosis, levels = c("CONTROL", "AS"))

# 5️) 논문용 박스플롯
p <- ggplot(pheno_data, aes(x = Diagnosis, y = OASL, fill = Diagnosis)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    color = "black",
    size = 0.8
  ) +
  geom_jitter(
    aes(color = Diagnosis),
    width = 0.2,
    size = 1.5,
    shape = 21,
    stroke = 0.8,
    fill = "white",
    alpha = 0.8
  ) +
  stat_compare_means(
    comparisons = list(c("CONTROL", "AS")),
    method = "wilcox.test",
    label = "p.signif",
    size = 5,
    tip.length = 0.01,
    label.y = max(pheno_data$OASL, na.rm = TRUE) * 1.05  # 자동 위치 조정
  ) +
  coord_cartesian(ylim = c(0, max(pheno_data$OASL, na.rm = TRUE) * 1.1)) +
  labs(
    title = "AS(GSE100927) - OASL Expression",
    x = "",
    y = "OASL Expression"
  ) +
  scale_fill_manual(values = c("CONTROL" = "#619CFF", "AS" = "#F8766D")) +
  scale_color_manual(values = c("CONTROL" = "#619CFF", "AS" = "#F8766D")) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# 플롯 출력
print(p)

# 저장하려면 이 코드도 추가
ggsave("Fig_OASL_Control_vs_AS_Boxplot.tiff", plot = p, device = "tiff", dpi = 600, width = 6, height = 5, units = "in")


























# 1️) mapped 결과에서 OASL에 해당하는 probe ID 찾기
oasl_probe <- mapped$agilent_sureprint_g3_ge_8x60k_v2[mapped$hgnc_symbol == "OASL"]

# 2️) 해당 probe의 발현값 추출 (log2 normalized matrix)
oasl_expr <- as.numeric(normalized_data[oasl_probe, ])

# 3️) OASL 발현 기준으로 그룹 나누기 (median 기준)
group <- ifelse(oasl_expr >= median(oasl_expr, na.rm = TRUE), "High", "Low")

# 4️) pheno_data에 그룹 정보 추가
pheno_data$OASL_Group <- group
table(pheno_data$OASL_Group)















# OASL_Group factor로 지정
pheno_data$OASL_Group <- factor(pheno_data$OASL_Group, levels = c("Low", "High"))

# 디자인 매트릭스 생성
design_oasl <- model.matrix(~0 + OASL_Group, data = pheno_data)
colnames(design_oasl) <- levels(pheno_data$OASL_Group)  # "Low", "High"

# 모델 피팅
fit_oasl <- lmFit(normalized_data, design_oasl)

# contrast 설정: High - Low
contrast_oasl <- makeContrasts(High - Low, levels = design_oasl)
fit_oasl2 <- contrasts.fit(fit_oasl, contrast_oasl)
fit_oasl2 <- eBayes(fit_oasl2)

# DEG 결과 추출
deg_oasl <- topTable(fit_oasl2, adjust = "BH", number = Inf)
deg_oasl$ID <- rownames(deg_oasl)  # probe ID 저장

# biomaRt로 gene symbol 매핑 (이미 `mapped`에 있음)
deg_oasl <- left_join(deg_oasl, mapped, by = c("ID" = "agilent_sureprint_g3_ge_8x60k_v2"))

# 결과 확인
head(deg_oasl)

# 저장 (선택)
write.csv(deg_oasl, "OASL_High_vs_Low_DEG_microarray.csv", row.names = FALSE)


# 분석 함수 정의
run_oasl_deg_by_group <- function(group_name, expr_matrix, pheno_data, oasl_probe) {
  
  # 1️) 그룹 내 샘플만 선택
  samples_in_group <- rownames(pheno_data)[pheno_data$Condition == group_name]
  expr_sub <- expr_matrix[, samples_in_group]
  pheno_sub <- pheno_data[samples_in_group, ]
  
  # 2️) OASL 발현 추출
  oasl_expr <- as.numeric(expr_sub[oasl_probe, ])
  
  # 3️) 발현 기준 그룹 나누기 (중앙값 기준)
  oasl_group <- ifelse(oasl_expr >= median(oasl_expr, na.rm = TRUE), "High", "Low")
  pheno_sub$OASL_Group <- factor(oasl_group, levels = c("Low", "High"))
  
  # 4️) 디자인 매트릭스 및 모델
  design <- model.matrix(~0 + OASL_Group, data = pheno_sub)
  colnames(design) <- levels(pheno_sub$OASL_Group)
  
  # 5️) 분석
  fit <- lmFit(expr_sub, design)
  contrast <- makeContrasts(High - Low, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  
  # 6️) 결과
  deg <- topTable(fit2, number = Inf, adjust = "BH")
  deg$ID <- rownames(deg)
  deg <- left_join(deg, mapped, by = c("ID" = "agilent_sureprint_g3_ge_8x60k_v2"))
  
  # 7️) 저장
  filename <- paste0("DEG_OASL_High_vs_Low_", group_name, ".csv")
  write.csv(deg, filename, row.names = FALSE)
  cat(group_name, ": 저장 완료 →", sum(deg$adj.P.Val < 0.05), "significant DEGs\n")
  
  return(deg)
}

# OASL probe ID (기존 mapped에서 추출)
oasl_probe <- mapped$agilent_sureprint_g3_ge_8x60k_v2[mapped$hgnc_symbol == "OASL"]

# 실행
deg_control <- run_oasl_deg_by_group("Control", normalized_data, pheno_data, oasl_probe)
deg_other   <- run_oasl_deg_by_group("Other", normalized_data, pheno_data, oasl_probe)


# 분석 함수 정의
run_oasl_deg_by_group <- function(group_name, expr_matrix, pheno_data, oasl_probe) {
  
  # 1️) 그룹 내 샘플만 선택
  samples_in_group <- rownames(pheno_data)[pheno_data$Condition == group_name]
  expr_sub <- expr_matrix[, samples_in_group]
  pheno_sub <- pheno_data[samples_in_group, ]
  
  # 2️) OASL 발현 추출
  oasl_expr <- as.numeric(expr_sub[oasl_probe, ])
  
  # 3️) 발현 기준 그룹 나누기 (중앙값 기준)
  oasl_group <- ifelse(oasl_expr >= median(oasl_expr, na.rm = TRUE), "High", "Low")
  pheno_sub$OASL_Group <- factor(oasl_group, levels = c("Low", "High"))
  
  # 4️) 디자인 매트릭스 및 모델
  design <- model.matrix(~0 + OASL_Group, data = pheno_sub)
  colnames(design) <- levels(pheno_sub$OASL_Group)
  
  # 5️) 분석
  fit <- lmFit(expr_sub, design)
  contrast <- makeContrasts(High - Low, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  
  # 6️) 결과
  deg <- topTable(fit2, number = Inf, adjust = "BH")
  deg$ID <- rownames(deg)
  deg <- left_join(deg, mapped, by = c("ID" = "agilent_sureprint_g3_ge_8x60k_v2"))
  
  # 7️) 저장
  filename <- paste0("DEG_OASL_High_vs_Low_", group_name, ".csv")
  write.csv(deg, filename, row.names = FALSE)
  cat(group_name, ": 저장 완료 →", sum(deg$adj.P.Val < 0.05), "significant DEGs\n")
  
  return(deg)
}

# OASL probe ID (기존 mapped에서 추출)
oasl_probe <- mapped$agilent_sureprint_g3_ge_8x60k_v2[mapped$hgnc_symbol == "OASL"]

# 실행
deg_control <- run_oasl_deg_by_group("Control", normalized_data, pheno_data, oasl_probe)
deg_other   <- run_oasl_deg_by_group("Other", normalized_data, pheno_data, oasl_probe)

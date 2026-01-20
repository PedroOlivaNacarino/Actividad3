#### ============================================================
####  
####  ACTIVIDAD 3
####  
#### ============================================================

#### Paquetes necesarios

rm(list = ls())
df<-read.csv("Base de datos Los Simpson completa.csv")
View(df)

#### Paquetes necesarios
library(tidyverse)
library(factoextra)
library(pheatmap)
library(gtsummary)
library(broom)
library(performance)
library(pROC)


#### ============================================================
####  PASO 1. Limpieza de datos
#### ============================================================

df_clean <- df %>% filter(complete.cases(select(., ADCY3:UHMK1)))
any(is.na(df_clean %>% select(ADCY3:UHMK1)))

View(df_clean)

df_genes <- df_clean%>% select(ADCY3:UHMK1)
View(df_genes)

colSums(is.na(df_genes))

#### ============================================================
####  PASO 2. Normalidad
#### ============================================================

shapiro_res <- df_genes %>%
  pivot_longer(everything(), names_to="gen", values_to="expr") %>%
  group_by(gen) %>%
  summarise(p_value = shapiro.test(expr)$p.value, .groups="drop") %>%
  mutate(normal_aprox = ifelse(p_value > 0.05, "SI", "NO"))

table(shapiro_res$normal_aprox)

norm_genes_testSha <- lapply(df_genes, shapiro.test)
norm_genes_testSha

resumen_normalidad <- norm_genes_testSha %>%
  map_df(~list(
    W_Statistic = .x$statistic, 
    P_Value = .x$p.value
  ), .id = "Gen")
resumen_normalidad

#### ============================================================
####  PASO 3. Análisis de Componentes Principales (PCA)
#### ============================================================

pca <- prcomp(df_genes, center = TRUE, scale. = TRUE)
summary(pca)

scores <- as.data.frame(pca$x[,1:6])
colnames(scores) <- paste0("PC",1:6)
df_pca <- bind_cols(df_clean, scores)
View(df_pca)


pca_res <- prcomp(
  df_genes %>% select(ADCY3:UHMK1),
  center = TRUE,
  scale. = TRUE
)

var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2) 

tabla_pca <- tibble(
  Componente = paste0("PC", 1:6),
  R2 = round(var_exp[1:6], 3)
)

tabla_pca_final <- tabla_pca %>%
  gt::gt() %>%
  gt::tab_caption(
    "Tabla X. Varianza explicada (R²) por las primeras seis componentes principales del PCA basado en la expresión génica"
  )

tabla_pca_final

#### ============================================================
####  PASO 4. Visualización y clustering
#### ============================================================

# Variables PCA coloreadas por cos2
fviz_pca_var(pca, col.var="cos2" ) +
  ggtitle("PCA - Variables (cos2)") +
  theme_classic()

# Distribucon de individuos en las dos primeras componentes principales del PCA (PC1 y PC2)
fviz_pca_ind(pca, geom="point") +
  ggtitle("PCA - Individuos") +
  theme_classic()

# Clusterizacion
km_pca <- kmeans(df_pca %>% select(PC1:PC6), centers = 3, nstart = 25)

#### Visualización del clustering
fviz_cluster(km_pca,
             data = df_pca %>% select(PC1:PC6),
             geom = "point",
             ellipse.type = "norm") +
  theme_classic()

# Contribucion de variables a las componentes
fviz_contrib(pca, choice = "var", axes = 1, top = 37, ggtheme = theme_minimal())
fviz_contrib(pca, choice = "var", axes = 2, top = 37, ggtheme = theme_minimal())

# Pacientes coloreados por IMC (PC1 - PC2)

df_pca <- df_pca %>%
  mutate(
    categoria_imc = case_when(
      imc_kg_m2 < 25 ~ "Normal",
      imc_kg_m2 >= 25 & imc_kg_m2 < 30 ~ "Sobrepeso",
      imc_kg_m2 >= 30 ~ "Obesidad"
    )
  )
df_pca


graficoIMC <- fviz_pca_ind(pca, geom.ind = "point",
                           col.ind = df_pca$categoria_imc,
                           palette = c("Normal" = "green", "Sobrepeso" = "orange", "Obesidad" = "red"),
                           addEllipses = TRUE,
                           legend.title = "Categoría IMC")+
  ggtitle("PCA - Individuos según categoría IMC") +
  theme_classic()

graficoIMC

# Clustering basado en scores de PCA (3 clusters)
score_pca <- pca$x[,1:2]

set.seed(2026)

km <- kmeans(score_pca, centers = 3, nstart = 50)

fviz_cluster(km, data = score_pca, geom = "point") +
  theme_classic()


#### ============================================================
####  PASO 5. Heatmap
#### ============================================================

cor_matriz <- cor(df_pca %>% select(ADCY3:UHMK1),
            df_pca %>% select(PC1:PC6),
            method="spearman")

cor_matriz <- cor_matriz[ apply(cor_matriz, 1, function(x) all(is.finite(x))), apply(cor_matriz, 2, function(x) all(is.finite(x)))]

pheatmap(cor_matriz, 
         clustering_method = "complete",
         fontsize_row = 9,
         fontsize_col= 10,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         angle_col = 45,
         border_color = "grey90",
         main="Correlación Spearman entre genes y componentes principales")

#### ============================================================
####  PASO 6. Tabla descriptiva
#### ============================================================

make_terciles <- function(x){
  q <- quantile(x, probs=c(1/3,2/3))
  cut(x, breaks=c(-Inf,q,Inf), labels=c("t1","t2","t3"), include.lowest = TRUE)
}

df_pca <- df_pca %>%
  mutate(
    PC1_t = make_terciles(PC1),
    PC2_t = make_terciles(PC2),
    PC3_t = make_terciles(PC3)
  )

tablaPC1 <- df_pca %>%
  select(PC1_t, imc_kg_m2, edad_anios, ADCY3:UHMK1) %>%
  tbl_summary(
    by = PC1_t,
    statistic = all_continuous() ~ "{median} ({p25}-{p75})",
    digits = all_continuous() ~ 2) %>%
  add_p(test = all_continuous() ~ "kruskal.test")

tablaPC2 <- df_pca %>%
  select(PC2_t, imc_kg_m2, edad_anios, ADCY3:UHMK1) %>%
  tbl_summary(
    by = PC2_t,
    statistic = all_continuous() ~ "{median} ({p25}-{p75})",
    digits = all_continuous() ~ 2) %>%
  add_p(test = all_continuous() ~ "kruskal.test")

tablaPC3 <- df_pca %>%
  select(PC3_t, imc_kg_m2, edad_anios, ADCY3:UHMK1) %>%
  tbl_summary(
    by = PC3_t,
    statistic = all_continuous() ~ "{median} ({p25}-{p75})",
    digits = all_continuous() ~ 2) %>%
  add_p(test = all_continuous() ~ "kruskal.test")

tablaPCfinal <- tbl_merge(
  tbls = list(tablaPC1, tablaPC2, tablaPC3),
  tab_spanner = c(
    "Terciles PC1",
    "Terciles PC2",
    "Terciles PC3"))

tablaPCfinal


#### ============================================================
####  PASO 7. Regresión logística
#### ============================================================

# Como no hay una columna que indique metastasis o no, se considern IL6, y TNF-alpha como candidatos a indicadores tumorales al estar involucrados en respuesta pro-inflamatoria.

# Con IL-6
p_il6 <- median(df_pca$il6_pg_ml, na.rm = TRUE)

df_pca <- df_pca %>%
  mutate(
    il6_cat = ifelse(il6_pg_ml < p_il6, "Baja", "Alta"),
    il6_cat = factor(il6_cat, levels = c("Baja", "Alta")))

table(df_pca$il6_cat)

df_pca <- df_pca %>%
  mutate(
    PC1_t = relevel(factor(PC1_t), ref = "t1"),
    PC2_t = relevel(factor(PC2_t), ref = "t1"),
    PC3_t = relevel(factor(PC3_t), ref = "t1"),
    sexo = factor(sexo))

modelo_il6_log <- glm(
  il6_cat ~ PC1_t + PC2_t + PC3_t + edad_anios + sexo + imc_kg_m2,
  data = df_pca,
  family = binomial)

summary(modelo_log)

install.packages("broom.helpers")
library(broom.helpers)

tabla_logisticaIL6 <- tbl_regression(
  modelo_il6_log,
  exponentiate = TRUE
) %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Regresión logística para IL-6 categorizado (OR y p-values)**")

tabla_logisticaIL6


# Con TNF-alfa
p_tnfalfa <- median(df_pca$tnf_alpha_pg_ml, na.rm = TRUE)

df_pca <- df_pca %>%
  mutate(
    tnf_cat = ifelse(tnf_alpha_pg_ml < p_tnfalfa, "Baja", "Alta"),
    tnf_cat = factor(tnf_cat, levels = c("Baja", "Alta")))

table(df_pca$tnf_cat)

df_pca <- df_pca %>%
  mutate(
    PC1_t = relevel(factor(PC1_t), ref = "t1"),
    PC2_t = relevel(factor(PC2_t), ref = "t1"),
    PC3_t = relevel(factor(PC3_t), ref = "t1"),
    sexo = factor(sexo))

modelo_tnf_log <- glm(
  tnf_cat ~ PC1_t + PC2_t + PC3_t + edad_anios + sexo + imc_kg_m2,
  data = df_pca,
  family = binomial)

summary(modelo_log)

tabla_logisticaTNF <- tbl_regression(
  modelo_tnf_log,
  exponentiate = TRUE
) %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Regresión logística para TNF-alpha categorizado (OR y p-values)**")

tabla_logisticaTNF


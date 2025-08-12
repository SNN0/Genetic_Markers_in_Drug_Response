
library(grid)
library(ggplot2)
library(gplots)
library(ggfortify)
library(ggpubr)
library(ggraph)
library(grid)
library(gridExtra)
library(survminer)
library(cowplot)

#km için çözünürlük : 715 X 310

TCF7L2_KM
TCF7L2_BRCA_KM

FLNB_KM
FLNB_BRCA_KM

LDHB_KM
LDHB_BRCA_KM

MMP7_KM
MMP7_BRCA_KM

FYN_KM
FYN_BRCA_KM

CXCR4_KM
CXCR6_BRCA_KM




TCF7L2_KM$plot$layers=TCF7L2_KM$plot$layers[-4]
TCF7L2_KM$plot$layers[[1]]$aes_params$size <- 0.8


# Her şeyi küçültmek için:
TCF7L2_KM$plot <- TCF7L2_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p < 0.0001",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel Samples with TCF7L2\n(CNA x Response, HR = 2.94, with 95% CI: [1, 5.73])',
               x='Days')

TCF7L2_KM


#-------

TCF7L2_BRCA_KM$plot$layers=TCF7L2_BRCA_KM$plot$layers[-4]
TCF7L2_BRCA_KM$plot$layers[[1]]$aes_params$size <- 0.8


# Her şeyi küçültmek için:
TCF7L2_BRCA_KM$plot <- TCF7L2_BRCA_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.83",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel BRCA Samples with TCF7L2\n(CNA x Response, HR = 1, with 95% CI: [0.99, 1])',
               x='Days')

TCF7L2_BRCA_KM


TCF7L2_ALLvsBRCA = grid.arrange(TCF7L2_KM$plot,TCF7L2_BRCA_KM$plot,ncol=2,nrow=1)



#------------
FLNB_KM



FLNB_KM$plot$layers=FLNB_KM$plot$layers[-4]
FLNB_KM$plot$layers[[1]]$aes_params$size <- 0.8

FLNB_KM$plot <- FLNB_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
               
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p < 0.0001",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel Samples with FLNB\n(CNA x Response, HR = 4.83, with 95% CI: [1, 19.54])',
               x='Days')

FLNB_KM


#------------
FLNB_BRCA_KM



FLNB_BRCA_KM$plot$layers=FLNB_BRCA_KM$plot$layers[-4]
FLNB_BRCA_KM$plot$layers[[1]]$aes_params$size <- 0.8

FLNB_BRCA_KM$plot <- FLNB_BRCA_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.45",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel BRCA Samples with FLNB\n(CNA x Response, HR = 1, with 95% CI: [1, 1])',
               x='Days')

FLNB_BRCA_KM



FLNB_ALLvsBRCA = grid.arrange(FLNB_KM$plot,FLNB_BRCA_KM$plot,ncol=2,nrow=1)



#LDHB_KM-------------------------------------------------
LDHB_KM

LDHB_KM$plot$layers=LDHB_KM$plot$layers[-4]
LDHB_KM$plot$layers[[1]]$aes_params$size <- 0.8

LDHB_KM$plot <- LDHB_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.00063",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel Samples with LDHB\n(CNA x Response, HR = 3.47, with 95% CI: [0.71, 12.96])',
               x='Days')



LDHB_KM



#LDHB_KM-------------------------------------------------
LDHB_BRCA_KM

LDHB_BRCA_KM$plot$layers=LDHB_BRCA_KM$plot$layers[-4]
LDHB_BRCA_KM$plot$layers[[1]]$aes_params$size <- 0.8

LDHB_BRCA_KM$plot <- LDHB_BRCA_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.83",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel BRCA Samples with LDHB\n(CNA x Response, HR = 1, with 95% CI: [0.99, 1])',
               x='Days')



LDHB_BRCA_KM

LDHB_ALLvsBRCA = grid.arrange(LDHB_KM$plot,LDHB_BRCA_KM$plot,ncol=2,nrow=1)



#MMP7_KM *-------------------------------------

MMP7_KM

MMP7_KM$plot$layers=MMP7_KM$plot$layers[-4]
MMP7_KM$plot$layers[[1]]$aes_params$size <- 0.8

MMP7_KM$plot <- MMP7_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.00016",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel Samples with MMP7\n(CNA x Response, HR = 3.05, with 95% CI: [1, 9.52])',
               x='Days')



MMP7_KM


#MMP7_KM *-------------------------------------

MMP7_BRCA_KM

MMP7_BRCA_KM$plot$layers=MMP7_BRCA_KM$plot$layers[-4]
MMP7_BRCA_KM$plot$layers[[1]]$aes_params$size <- 0.8

MMP7_BRCA_KM$plot <- MMP7_BRCA_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.66",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel Samples with MMP7\n(CNA x Response, HR = 0.99, with 95% CI: [0.98, 1])',
               x='Days')



MMP7_BRCA_KM




MMP7_ALLvsBRCA = grid.arrange(MMP7_KM$plot,MMP7_BRCA_KM$plot,ncol=2,nrow=1)



#FYN_KM *-------------------------------------

FYN_KM

FYN_KM$plot$layers=FYN_KM$plot$layers[-4]
FYN_KM$plot$layers[[1]]$aes_params$size <- 0.8

FYN_KM$plot <- FYN_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.00094",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel Samples with FYN\n(CNA x Response, HR = 0.21, with 95% CI: [0.03, 2.31])',
               x='Days')



FYN_KM


#MMP7_KM *-------------------------------------

FYN_BRCA_KM

FYN_BRCA_KM$plot$layers=FYN_BRCA_KM$plot$layers[-4]
FYN_BRCA_KM$plot$layers[[1]]$aes_params$size <- 0.8

FYN_BRCA_KM$plot <- FYN_BRCA_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
               
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.72",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel BRCA Samples with FYN\n(CNA x Response, HR = 1, with 95% CI: [0.99, 1])',
               x='Days')



FYN_BRCA_KM




FYN_ALLvsBRCA = grid.arrange(FYN_KM$plot,FYN_BRCA_KM$plot,ncol=2,nrow=1)


#CXCR4_KM *-------------------------------------

CXCR4_KM

CXCR4_KM$plot$layers=CXCR4_KM$plot$layers[-4]
CXCR4_KM$plot$layers[[1]]$aes_params$size <- 0.8

CXCR4_KM$plot <- CXCR4_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.004",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel Samples with CXCR4\n(CNA x Response, HR = 1, with 95% CI: [0.99, 1])',
               x='Days')



CXCR4_KM


#CXCR4_BRCA_KM *-------------------------------------

CXCR4_BRCA_KM

CXCR4_BRCA_KM$plot$layers=CXCR4_BRCA_KM$plot$layers[-4]
CXCR4_BRCA_KM$plot$layers[[1]]$aes_params$size <- 0.8

CXCR4_BRCA_KM$plot <- CXCR4_BRCA_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.26",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Paclitaxel BRCA Samples with CXCR4\n(CNA x Response, HR = 1, with 95% CI: [1, 1])',
               x='Days')



CXCR4_BRCA_KM




CXCR4_ALLvsBRCA = grid.arrange(CXCR4_KM$plot,CXCR4_BRCA_KM$plot,ncol=2,nrow=1)


#--------------------------

CXCR4_KM
FYN_KM
MMP7_KM
LDHB_KM
FLNB_KM
TCF7L2_KM

Figure4 <- plot_grid(
        TCF7L2_KM$plot, 
        FLNB_KM$plot,
        LDHB_KM$plot,
        MMP7_KM$plot,
        FYN_KM$plot,
        CXCR4_KM$plot,
        ncol = 3,
        nrow=2,
        labels = c("a", "b",'c','d','e','f') # Etiketleri 
)


#-----------------------------

#BAR PLOT çözünürlük:419 X 576 

CXCR4_FYN_HR_BRCA

CXCR4_FYN_HR_BRCA$labels$title = element_blank()
CXCR4_FYN_HR_BRCA$theme$axis.title.y$size = 6
CXCR4_FYN_HR_BRCA$theme$axis.title.y$face='bold'
CXCR4_FYN_HR_BRCA$theme$strip.text$size=7
CXCR4_FYN_HR_BRCA$theme$strip.text$face='bold'
CXCR4_FYN_HR_BRCA$theme$axis.text.x$size=5
CXCR4_FYN_HR_BRCA$theme$axis.text.y$size=6
CXCR4_FYN_HR_BRCA$theme$plot.margin = unit(c(-0.1, 0.1, 0.1, 0.1), "cm")


plot <- CXCR4_FYN_HR_BRCA
plot$layers[[2]]$aes_params$size <- 1.8  # Yazı boyutunu küçültmek için 3 olarak 
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak 


CXCR4_FYN_HR_BRCA = plot
CXCR4_FYN_HR_BRCA

#C index
CXCR4_FYN_BRCA_Cindex


CXCR4_FYN_BRCA_Cindex$labels$title = element_blank()
CXCR4_FYN_BRCA_Cindex$theme$axis.title.y$size = 7
CXCR4_FYN_BRCA_Cindex$theme$axis.text.y$size = 6
CXCR4_FYN_BRCA_Cindex$theme$axis.text.x$size= 7
CXCR4_FYN_BRCA_Cindex$theme$axis.title.y$face='bold'
CXCR4_FYN_BRCA_Cindex$theme$strip.text$size = 8
CXCR4_FYN_BRCA_Cindex$theme$plot.margin = unit(c(-0.1, 0.1, 0.1, 0.1), "cm")

plot <- CXCR4_FYN_BRCA_Cindex
plot$layers[[1]]$geom_params$width <- 0.5
plot$layers[[2]]$aes_params$size<- 2.4# Yazı boyutunu küçültmek için 3 olarak 
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak 
#print(plot)

CXCR4_FYN_BRCA_Cindex = plot

CXCR4_FYN_BRCA_Cindex


title <- textGrob(
        "Paclitaxel BRCA Samples (HR and C-Index of FYN & CXCR4)", 
        gp = gpar(fontsize = 9, fontface = "bold")  # Yazı boyutu ve stili
)

#1. grafik
TCF7L2_LDHB_Cindex_HR_Graph=grid.arrange(CXCR4_FYN_HR_BRCA,CXCR4_FYN_BRCA_Cindex,nrow=2,top=title)

plot(TCF7L2_LDHB_Cindex_HR_Graph)





#2. grafik # 700 x 330
FYN_TCF7L2_BRCA_Cindex_path

FYN_TCF7L2_BRCA_Cindex_path$labels$title = 'Paclitaxel BRCA Samples (C-Index)'
FYN_TCF7L2_BRCA_Cindex_path$theme$plot.title$size=14
FYN_TCF7L2_BRCA_Cindex_path$theme$text$face='bold'
FYN_TCF7L2_BRCA_Cindex_path$theme$axis.title.y$size = 9
FYN_TCF7L2_BRCA_Cindex_path$theme$axis.text.y$size = 9
FYN_TCF7L2_BRCA_Cindex_path$theme$axis.text.x$size= 9
FYN_TCF7L2_BRCA_Cindex_path$theme$axis.title.y$face='bold'
FYN_TCF7L2_BRCA_Cindex_path$theme$strip.text$size = 9

plot <- FYN_TCF7L2_BRCA_Cindex_path
plot$layers[[1]]$geom_params$width <- 0.7
plot$layers[[2]]$aes_params$size<- 3# Yazı boyutunu küçültmek için 3 olarak 
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak 
#print(plot)

FYN_TCF7L2_BRCA_Cindex_path = plot
FYN_TCF7L2_BRCA_Cindex_path



FYN_MYL9_FLNB_ALL_Cindex
FYN_MYL9_FLNB_BRCA_Cindex

FYN_TCF7L2_ALL_Cindex_path
FYN_TCF7L2_BRCA_Cindex_path


Figure2 <- plot_grid(
        FYN_MYL9_FLNB_ALL_Cindex, 
        FYN_MYL9_FLNB_BRCA_Cindex,
        FYN_TCF7L2_ALL_Cindex_path,
        FYN_TCF7L2_BRCA_Cindex_path,
        ncol = 2,
        nrow=2,
        labels = c("a", "b",'c','d') # Etiketleri 
)




#barlar



#4. grafik # 450 x 570
LYN_Cindex

LYN_Cindex$labels$title = element_blank()
LYN_Cindex$theme$axis.title.y$size = 9
LYN_Cindex$theme$axis.title.y$face='bold'
LYN_Cindex$theme$strip.text$size=8
LYN_Cindex$theme$strip.text$face='bold'
LYN_Cindex$theme$axis.text.x$size=7
LYN_Cindex$theme$axis.text.y$size=7
LYN_Cindex$theme$plot.margin = unit(c(-0.1, 0.1, 0.1, 0.1), "cm")

plot <- LYN_Cindex
#plot$layers[[1]]$geom_params$width <- 0.5
plot$layers[[2]]$aes_params$size<- 2.3# Yazı boyutunu küçültmek için 3 olarak 
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak 
#print(plot)

LYN_Cindex=plot

TNFSF13B_Cindex

TNFSF13B_Cindex$labels$title = element_blank()
TNFSF13B_Cindex$theme$axis.title.y$size = 9
TNFSF13B_Cindex$theme$axis.title.y$face='bold'
TNFSF13B_Cindex$theme$strip.text = element_blank()
TNFSF13B_Cindex$theme$strip.text$face='bold'
TNFSF13B_Cindex$theme$axis.text.x$size=7
TNFSF13B_Cindex$theme$axis.text.y$size=7
TNFSF13B_Cindex$theme$plot.margin = unit(c(-0.1, 0.1, 0.1, 0.1), "cm")


plot <- TNFSF13B_Cindex
#plot$layers[[1]]$geom_params$width <- 0.5
plot$layers[[2]]$aes_params$size<- 2.3# Yazı boyutunu küçültmek için 3 olarak 
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak
#print(plot)
TNFSF13B_Cindex=plot
TNFSF13B_Cindex


title2 <- textGrob(
        "5-FU STAD Samples of LYN and TNFSF13B (C-Index)", 
        gp = gpar(fontsize = 11, fontface = "bold")  # Yazı boyutu ve stili
)



TNFSF13B_Cindex_AND_LYN_Cindex=grid.arrange(LYN_Cindex,TNFSF13B_Cindex, top = title2)

plot(TNFSF13B_Cindex_AND_LYN_Cindex)





LYN_Cindex_ALL
TNFSF13B_Cindex_ALL

#450 x 570

LYN_Cindex_ALL

LYN_Cindex_ALL$labels$title = element_blank()
LYN_Cindex_ALL$theme$axis.title.y$size = 9
LYN_Cindex_ALL$theme$axis.title.y$face='bold'
LYN_Cindex_ALL$theme$strip.text$size=8
LYN_Cindex_ALL$theme$strip.text$face='bold'
LYN_Cindex_ALL$theme$axis.text.x$size=7
LYN_Cindex_ALL$theme$axis.text.y$size=7
LYN_Cindex_ALL$theme$plot.margin = unit(c(-0.1, 0.1, 0.1, 0.1), "cm")

plot <- LYN_Cindex_ALL
#plot$layers[[1]]$geom_params$width <- 0.5
plot$layers[[2]]$aes_params$size<- 2.3# Yazı boyutunu küçültmek için 3 olarak 
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak 
#print(plot)

LYN_Cindex_ALL=plot

TNFSF13B_Cindex_ALL

TNFSF13B_Cindex_ALL$labels$title = element_blank()
TNFSF13B_Cindex_ALL$theme$axis.title.y$size = 9
TNFSF13B_Cindex_ALL$theme$axis.title.y$face='bold'
TNFSF13B_Cindex_ALL$theme$strip.text = element_blank()
TNFSF13B_Cindex_ALL$theme$strip.text$face='bold'
TNFSF13B_Cindex_ALL$theme$axis.text.x$size=7
TNFSF13B_Cindex_ALL$theme$axis.text.y$size=7
TNFSF13B_Cindex_ALL$theme$plot.margin = unit(c(-0.1, 0.1, 0.1, 0.1), "cm")


plot <- TNFSF13B_Cindex_ALL
#plot$layers[[1]]$geom_params$width <- 0.5
plot$layers[[2]]$aes_params$size<- 2.3# Yazı boyutunu küçültmek için 3 olarak 
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak 
#print(plot)
TNFSF13B_Cindex_ALL=plot
TNFSF13B_Cindex_ALL


title2 <- textGrob(
        "5-FU Samples of LYN and TNFSF13B (C-Index)", 
        gp = gpar(fontsize = 11, fontface = "bold")  # Yazı boyutu ve stili
)



TNFSF13B_Cindex_AND_LYN_Cindex_ALL=grid.arrange(LYN_Cindex_ALL,TNFSF13B_Cindex_ALL, top = title2)

plot(TNFSF13B_Cindex_AND_LYN_Cindex_ALL)









#5. grafik   870 x 360
LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex=LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex$CNA

LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex$labels$title = '5-FU STAD Samples (C-Index)'

LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex$theme$plot.title$size=12
LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex$theme$axis.title.y$size = 9
LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex$theme$axis.title.y$face='bold'
LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex$theme$strip.text$size=10
LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex$theme$strip.text$face='bold'
LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex$theme$axis.text.x$size=7
LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex$theme$axis.text.y$size=8
LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex$theme$plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")

plot <- LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex
#plot$layers[[1]]$geom_params$width <- 0.5
plot$layers[[2]]$aes_params$size<- 2.3# Yazı boyutunu küçültmek için 3 olara
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak 
#print(plot)
LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex=plot
LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex



# ALL 5FU
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex = LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex$CNA


LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex$labels$title = '5-FU Samples (C-Index)'

LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex$theme$plot.title$size=12
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex$theme$axis.title.y$size = 9
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex$theme$axis.title.y$face='bold'
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex$theme$strip.text$size=10
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex$theme$strip.text$face='bold'
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex$theme$axis.text.x$size=7
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex$theme$axis.text.y$size=8
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex$theme$plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")

plot <- LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex
#plot$layers[[1]]$geom_params$width <- 0.5
plot$layers[[2]]$aes_params$size<- 2.3# Yazı boyutunu küçültmek için 3 olara
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak 
#print(plot)
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex=plot
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex



LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex
LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex


plot(TNFSF13B_Cindex_AND_LYN_Cindex)
plot(TNFSF13B_Cindex_AND_LYN_Cindex_ALL)

plot(LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex)

Figure5 <- plot_grid(
        LYN_MYO1F_PINK1_TNFSF13B_PREX1_ALL_Cindex, 
        LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex,
        TNFSF13B_Cindex_AND_LYN_Cindex_ALL,
        TNFSF13B_Cindex_AND_LYN_Cindex,
        ncol = 2,
        nrow=2,
        labels = c("A", "B",'C','D'), # Etiketleri 
        rel_heights = c(1, 2),
        align = 'v' # Dikey olarak hizala
)

Figure5 <- plot_grid(
        TNFSF13B_Cindex_AND_LYN_Cindex_ALL,
        TNFSF13B_Cindex_AND_LYN_Cindex,
        ncol = 2,
        nrow=1,
        labels = c("C", "D") # Etiketleri 
        
        
)





#------------






#5FU - KM   715 X 310


BAMBI_KM_ALLSamples



BAMBI_KM_ALLSamples$plot$layers=BAMBI_KM_ALLSamples$plot$layers[-4]
BAMBI_KM_ALLSamples$plot$layers[[1]]$aes_params$size <- 0.8

BAMBI_KM_ALLSamples$plot <- BAMBI_KM_ALLSamples$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 1),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 1)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.00027",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='5-FU Samples with BAMBI\n(Mutation, HR = 7.28, with 95% CI: [1, 13.02])',
               x='Days')

BAMBI_KM_ALLSamples


#--------------------------------------


CDH1_5FU_ALL_KM


CDH1_5FU_ALL_KM$plot$layers=CDH1_5FU_ALL_KM$plot$layers[-4]
CDH1_5FU_ALL_KM$plot$layers[[1]]$aes_params$size <- 0.8

CDH1_5FU_ALL_KM$plot <- CDH1_5FU_ALL_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 1),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 1)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.39",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='5-FU Samples with CDH1\n(Mutation, HR = 1.01, with 95% CI: [1, 1])',
               x='Days')

CDH1_5FU_ALL_KM

#----

CDH1_KM_STADSamples


CDH1_KM_STADSamples$plot$layers=CDH1_KM_STADSamples$plot$layers[-4]
CDH1_KM_STADSamples$plot$layers[[1]]$aes_params$size <- 0.8

CDH1_KM_STADSamples$plot <- CDH1_KM_STADSamples$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 1),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 1)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.044",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='5-FU STAD Samples with CDH1\n(Mutation, HR = 1.18, with 95% CI: [1, 1.22])',
               x='Days')

CDH1_KM_STADSamples



CDH1_ALLvsSTAD =grid.arrange(CDH1_5FU_ALL_KM$plot,CDH1_KM_STADSamples$plot,ncol=2)




#



#----

LIMA1_5FU_ALL_KM


LIMA1_5FU_ALL_KM$plot$layers=LIMA1_5FU_ALL_KM$plot$layers[-4]
LIMA1_5FU_ALL_KM$plot$layers[[1]]$aes_params$size <- 0.8

LIMA1_5FU_ALL_KM$plot <- LIMA1_5FU_ALL_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.026",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='5-FU Samples with LIMA1\n(CNA x Response, HR = 1.98, with 95% CI: [0.89, 6.93])',
               x='Days')

LIMA1_5FU_ALL_KM







#*---

SIX4_STAD_KM


SIX4_STAD_KM$plot$layers=SIX4_STAD_KM$plot$layers[-4]
SIX4_STAD_KM$plot$layers[[1]]$aes_params$size <- 0.8

SIX4_STAD_KM$plot <- SIX4_STAD_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 1),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 1)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.0039",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='5-FU STAD Samples with SIX4\n(CNA, HR = 2.68, with 95% CI: [1, 5.16])',
               x='Days')

SIX4_STAD_KM



LIMA1_ALLvsSTAD =grid.arrange(LIMA1_5FU_ALL_KM$plot,LIMA1_5FU_STAD_KM$plot,ncol=2)



#715 X 310





KRT80_5FU_ALL_KM


#----

KRT80_5FU_ALL_KM


KRT80_5FU_ALL_KM$plot$layers=KRT80_5FU_ALL_KM$plot$layers[-4]
KRT80_5FU_ALL_KM$plot$layers[[1]]$aes_params$size <- 0.8

KRT80_5FU_ALL_KM$plot <- KRT80_5FU_ALL_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.043",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='5-FU Samples with KRT80\n(CNA x Response, HR = 1.21, with 95% CI: [1, 1.26])',
               x='Days')

KRT80_5FU_ALL_KM



#------------
FRMD6_STAD_KM

FRMD6_STAD_KM$plot$layers=FRMD6_STAD_KM$plot$layers[-4]
FRMD6_STAD_KM$plot$layers[[1]]$aes_params$size <- 0.8

FRMD6_STAD_KM$plot <- FRMD6_STAD_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 1),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 1)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.0039",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='5-FU STAD Samples with FRMD6\n(CNA, HR = 5.93, with 95% CI: [1, 16.1])',
               x='Days')

FRMD6_STAD_KM


#715 310
KRT80_ALLvsSTAD =grid.arrange(KRT80_5FU_ALL_KM$plot,KRT80_5FU_STAD_KM$plot,ncol=2)


#--------------------




Figure6 <- plot_grid(
        BAMBI_KM_ALLSamples$plot, 
        LIMA1_5FU_ALL_KM$plot,
        KRT80_5FU_ALL_KM$plot,
        CDH1_KM_STADSamples$plot,
        SIX4_STAD_KM$plot,
        FRMD6_STAD_KM$plot,
        ncol = 3,
        nrow=2,
        labels = c("a", "b",'c','d','e','f') # Etiketleri 
        #rel_heights = c(1, 2),
        #align = 'v' # Dikey olarak hizala
)












#-------------





#gemcitabine #km için çözünürlük : 679 X 419

#6. grafil

MGST2_CDH1_BMP4_GPX8_Cindex

MGST2_CDH1_BMP4_GPX8_Cindex$labels$title='Gemcitabine PAAD Samples (C-Index)'
MGST2_CDH1_BMP4_GPX8_Cindex$theme$plot.title$size=10
MGST2_CDH1_BMP4_GPX8_Cindex$theme$axis.title.y$size = 6
MGST2_CDH1_BMP4_GPX8_Cindex$theme$axis.title.y$face='bold'
MGST2_CDH1_BMP4_GPX8_Cindex$theme$strip.text$size=8
MGST2_CDH1_BMP4_GPX8_Cindex$theme$strip.text$face='bold'
MGST2_CDH1_BMP4_GPX8_Cindex$theme$axis.text.x$size=5
MGST2_CDH1_BMP4_GPX8_Cindex$theme$axis.text.y$size=6
MGST2_CDH1_BMP4_GPX8_Cindex$theme$plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")

plot <- MGST2_CDH1_BMP4_GPX8_Cindex
plot$layers[[2]]$aes_params$size <- 2  # Yazı boyutunu küçültmek için 3 olarak
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak 
#print(plot)


MGST2_CDH1_BMP4_GPX8_Cindex = plot




MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL

MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL$labels$title='Gemcitabine Samples (C-Index)'
MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL$theme$plot.title$size=10
MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL$theme$axis.title.y$size = 6
MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL$theme$axis.title.y$face='bold'
MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL$theme$strip.text$size=8
MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL$theme$strip.text$face='bold'
MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL$theme$axis.text.x$size=5
MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL$theme$axis.text.y$size=6
MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL$theme$plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")

plot <- MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL
plot$layers[[2]]$aes_params$size <- 2  # Yazı boyutunu küçültmek için 3 olarak
plot$layers[[2]]$aes_params$vjust<- -0.2  # Yazı boyutunu küçültmek için 3 olarak 
#print(plot)

#680 450 x  
MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL = plot
library(cowplot)

MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL_ALLvsSTAD =grid.arrange(MGST2_CDH1_BMP4_GPX8_Cindex,MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL,ncol=1,nrow=2)

MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL_ALLvsSTAD <- plot_grid(
        MGST2_CDH1_BMP4_GPX8_Cindex_gemcitabine_ALL, 
        MGST2_CDH1_BMP4_GPX8_Cindex, 
        ncol = 1, 
        labels = c("A", "B") # Etiketleri 
)












#7. grafik

TBC1D10C_KM_PAAD

TBC1D10C_KM_PAAD$plot$layers=TBC1D10C_KM_PAAD$plot$layers[-4]
TBC1D10C_KM_PAAD$plot$layers[[1]]$aes_params$size <- 0.8

TBC1D10C_KM_PAAD$plot <- TBC1D10C_KM_PAAD$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.00061",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Gemcitabine PAAD Samples with TBC1D10C\n(Mutation x Response, HR = 2.62 with 95% CI: [1, 3.17])',
               x='Days')

TBC1D10C_KM_PAAD



#TBC1D10C_GEM_ALL_KM

TBC1D10C_GEM_ALL_KM

TBC1D10C_GEM_ALL_KM$plot$layers=TBC1D10C_GEM_ALL_KM$plot$layers[-4]
TBC1D10C_GEM_ALL_KM$plot$layers[[1]]$aes_params$size <- 0.8

TBC1D10C_GEM_ALL_KM$plot <- TBC1D10C_GEM_ALL_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.86",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Gemcitabine Samples with TBC1D10C\n(Mutation x Response, HR = 1.14 with 95% CI: [0.93, 1.65])',
               x='Days')

TBC1D10C_GEM_ALL_KM



TBC1D10C_ALLvsPAAD =grid.arrange(TBC1D10C_GEM_ALL_KM$plot,TBC1D10C_KM_PAAD$plot,ncol=2)








#LSR_KM_PAAD
LSR_KM_PAAD


LSR_KM_PAAD$plot$layers=LSR_KM_PAAD$plot$layers[-4]
LSR_KM_PAAD$plot$layers[[1]]$aes_params$size <- 0.8

LSR_KM_PAAD$plot <- LSR_KM_PAAD$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 2e-04",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Gemcitabine PAAD Samples with LSR\n(CNA x Response, HR = 4.15 with 95% CI: [1, 8.53])',
               x='Days')

LSR_KM_PAAD



#LSR_GEM_ALL_KM

LSR_GEM_ALL_KM

LSR_GEM_ALL_KM$plot$layers=LSR_GEM_ALL_KM$plot$layers[-4]
LSR_GEM_ALL_KM$plot$layers[[1]]$aes_params$size <- 0.8

LSR_GEM_ALL_KM$plot <- LSR_GEM_ALL_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.46",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Gemcitabine Samples with LSR\n(CNA x Response, HR = 1 with 95% CI: [0.99, 1])',
               x='Days')

LSR_GEM_ALL_KM


LSR_ALLvsPAAD =grid.arrange(LSR_GEM_ALL_KM$plot,LSR_KM_PAAD$plot,ncol=2)



#NQO1_KM_ALL

NQO1_KM_ALL


NQO1_KM_ALL$plot$layers=NQO1_KM_ALL$plot$layers[-4]
NQO1_KM_ALL$plot$layers[[1]]$aes_params$size <- 0.8

NQO1_KM_ALL$plot <- NQO1_KM_ALL$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p < 0.0001",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Gemcitabine Samples with NQO1\n(Mutation x Response, HR = 1.3 with 95% CI: [1, 1.30])',
               x='Days')

NQO1_KM_ALL




#NQO1_Gem_PAAD_KM
NQO1_Gem_PAAD_KM


NQO1_Gem_PAAD_KM$plot$layers=NQO1_Gem_PAAD_KM$plot$layers[-4]
NQO1_Gem_PAAD_KM$plot$layers[[1]]$aes_params$size <- 0.8

NQO1_Gem_PAAD_KM$plot <- NQO1_Gem_PAAD_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow = 2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.017",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Gemcitabine PAAD Samples with NQO1\n(Mutation x Response, HR = 1 with 95% CI: [0.99, 1])',
               x='Days')

NQO1_Gem_PAAD_KM



NQO1_ALLvsPAAD =grid.arrange(NQO1_KM_ALL$plot,NQO1_Gem_PAAD_KM$plot,ncol=2)







#CLIC3_KM_ALL

CLIC3_KM_ALL



CLIC3_KM_ALL$plot$layers=CLIC3_KM_ALL$plot$layers[-4]
CLIC3_KM_ALL$plot$layers[[1]]$aes_params$size <- 0.8

CLIC3_KM_ALL$plot <- CLIC3_KM_ALL$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow =2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p < 0.0001",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Gemcitabine Samples with CLIC3\n(CNA x Response, HR = 8.66 with 95% CI: [2.89, 29.20])',
               x='Days')

CLIC3_KM_ALL



#CLIC3_Gem_PAAD_KM

CLIC3_Gem_PAAD_KM



CLIC3_Gem_PAAD_KM$plot$layers=CLIC3_Gem_PAAD_KM$plot$layers[-4]
CLIC3_Gem_PAAD_KM$plot$layers[[1]]$aes_params$size <- 0.8

CLIC3_Gem_PAAD_KM$plot <- CLIC3_Gem_PAAD_KM$plot +
        theme_survminer() +
        theme(
                plot.title = element_text(size = 10, face = "bold",hjust = 0.5),  # Başlık küçültme
                axis.title.y = element_text(size = 9, face='bold'),  # Eksen başlıkları küçültme
                axis.title.x = element_text(size = 9,face='bold'),  # Eksen başlıkları küçültme
                axis.text.y = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme
                axis.text.x = element_text(size = 8,face='bold'),  # Eksen değerlerini küçültme),  # Eksen değerlerini küçültme
                legend.title = element_blank(),  # Legend başlığı
                legend.text = element_text(size = 7),  # Legend içeriği
                tex = element_text(size = 8),   # F
                legend.position = "top",
                # 
        ) +
        guides(
                color = guide_legend(ncol = 2, nrow =2),  # Renkler için sütun ve satır ayarı
                fill = guide_legend(ncol = 2, nrow = 2)   # Dolgular için sütun ve satır ayarı
        )+annotate(
                "text", 
                x = 500, y = 0.2,          # P-value yazısının pozisyonu
                label = "p = 0.47",      # Mevcut p-value
                size = 3,
                fontface='bold'# Font boyutu
        )+labs(title='Gemcitabine PAAD Samples with CLIC3\n(CNA x Response, HR = 1.09 with 95% CI: [0.88, 1.27])',
               x='Days')

CLIC3_Gem_PAAD_KM


TBC1D10C_ALLvsPAAD
LSR_ALLvsPAAD
NQO1_ALLvsPAAD
CLIC3_ALLvsPAAD 

CLIC3_ALLvsPAAD =grid.arrange(CLIC3_KM_ALL$plot,CLIC3_Gem_PAAD_KM$plot,ncol=2)


Figure7 <- plot_grid(
        NQO1_KM_ALL$plot, 
        CLIC3_KM_ALL$plot,
        TBC1D10C_KM_PAAD$plot,
        LSR_KM_PAAD$plot ,
        ncol = 2,
        nrow=2,
        labels = c("a", "b",'c','d') # Etiketleri 
)







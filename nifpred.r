    #setwd("C:/xampp/htdocs/nifPred/rcode/github")
    library(Biostrings)
    library(BioSeqClass)
    library(protr)
    library(e1071)
    library(R2HTML)
	
    if (file.exists("nifpred.html")) file.remove("nifpred.html")
	if (file.exists("nifpred.txt")) file.remove("nifpred.txt")
	x <- readAAStringSet("test.fasta")
    xx <- toupper(as.character(as.character(x)))
	###########Sequence filtering##########
	#Checking of standard residues#
std <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
zx <- xx
z <-  sapply(zx, function(s) strsplit(s, split=""))
zz <- lapply(z, table)
zz1 <- unique(as.character(unlist(lapply(zz, names))))
if(length(union(zz1, std))>20){
stop(" Could not execute the code because the sequence contains non-standard residues or ambigious character,so kindly submit sequences having standard residues only")
}else{
	###########End#################
	
    zm_ctd <- featureCTD(xx, class = elements("aminoacid"))
	res <- matrix(0, nrow = length(x), ncol = 2)
    tst <- data.frame(zm_ctd)
	s <- rep("u", nrow(tst))
    ts <- cbind(s, tst)
    SVM_rbf <- readRDS("svm_binary.rds")
    kk <- predict(SVM_rbf, newdata = ts[, -1], probability = TRUE)
    kz <- attr(kk, "probabilities")
    id <- which(kz[, "Y"] > 0.5)
    idn <- which(kz[, "Y"] <= 0.5)
    res[idn, 2] <- round(kz[, "N"][idn], 3)
    res[idn, 1] <- rep("non-nif", length(idn))
	
	mm <- which(!res[, 1] == "non-nif")
	
   
   if (length(mm) == 0) {
        res_f <- res
		resff <- data.frame(names(x),res_f)
	    colnames(resff) <- c("Sequence identifier", "Predicted as", "with probability")
	    HTML(resff, file = "nifpred.html", align = "center", Border = 1, innerBorder = 1)
		write.table(resff, "nifpred.txt", quote=F, sep="\t")
    }else {
        ts1 <- ts[id, ]
        SVM_rbf <- readRDS("jacknife.rds")
        kk <- predict(SVM_rbf, newdata = ts1[, -1])
        kz <- predict(SVM_rbf, newdata = ts1[, -1], probability = TRUE)
        zz <- round(attr(kz, "probabilities"), 3)
        res[id, 1] <- as.character(kk)
        res[id, 2] <- apply(zz, 1, max)
        res_f <- res
    	resff <- data.frame(names(x),res_f)
		colnames(resff) <- c("Sequence identifier", "Predicted as", "with probability")
		write.table(resff, "nifpred.txt", quote=F, sep="\t")
        kk <- which(!res_f[, 1] == "non-nif")
        kx <- x[kk]
        zz <- res[kk, ]
		if (length(kk) == 1){
		res_final <- data.frame(names(kx), zz[1], zz[2])
		colnames(res_final) <- c("Sequence identifier", "Predicted as", "with probability")
	    HTML(res_final, file = "nifpred.html", align = "center", Border = 1, innerBorder = 1)
		
		}else{
        zf <- zz[order(data.frame(zz[, 1])), ]
        zx <- kx[order(data.frame(zz[, 1]))]
        mf <- which(as.numeric(zf[, 2]) >= 0.4)
        resf <- zf[mf, ]
        seqf <- names(zx[mf])
        res_ff <- data.frame(seqf, resf)
		
		
        if (length(which(res_ff[, 2] == "nifB")) > 3) {
            mp_b <- res_ff[which(res_ff[, 2] == "nifB"), ]
            mpb <- mp_b[order(data.frame(mp_b[, 3]), decreasing = TRUE),]
            mzb <- mpb[1:3, ]
			
        }else {
            mzb <- res_ff[which(res_ff[, 2] == "nifB"), ]
        }
        
		if (length(which(res_ff[, 2] == "nifD")) > 3) {
            mp_d <- res_ff[which(res_ff[, 2] == "nifD"), ]
            mpd <- mp_d[order(data.frame(mp_d[, 3]), decreasing = TRUE),]
            mzd <- mpd[1:3, ]
        }else {
            mzd <- res_ff[which(res_ff[, 2] == "nifD"), ]
        }
        if (length(which(res_ff[, 2] == "nifE")) > 3) {
            mp_e <- res_ff[which(res_ff[, 2] == "nifE"), ]
            mpe <- mp_e[order(data.frame(mp_e[, 3]), decreasing = TRUE),]
            mze <- mpe[1:3, ]
        }else {
            mze <- res_ff[which(res_ff[, 2] == "nifE"), ]
        }
        
		if (length(which(res_ff[, 2] == "nifH")) > 3) {
            mp_h <- res_ff[which(res_ff[, 2] == "nifH"), ]
            mph <- mp_h[order(data.frame(mp_h[, 3]), decreasing = TRUE),]
            mzh <- mph[1:3, ]
        }else {
            mzh <- res_ff[which(res_ff[, 2] == "nifH"), ]
        }
       
	   if (length(which(res_ff[, 2] == "nifK")) > 3) {
            mp_k <- res_ff[which(res_ff[, 2] == "nifK"), ]
            mpk <- mp_k[order(data.frame(mp_k[, 3]), decreasing = TRUE),]
            mzk <- mpk[1:3, ]
        }else {
            mzk <- res_ff[which(res_ff[, 2] == "nifK"), ]
        }
        
		if (length(which(res_ff[, 2] == "nifN")) > 3) {
            mp_n <- res_ff[which(res_ff[, 2] == "nifN"), ]
            mpn <- mp_n[order(data.frame(mp_n[, 3]), decreasing = TRUE), ]
            mzn <- mpn[1:3, ]
       }else {
            mzn <- res_ff[which(res_ff[, 2] == "nifN"), ]
        }
        res_final <- data.frame(rbind(mzb, mzd, mze, mzh, mzk, mzn))
		colnames(res_final) <- c("Sequence identifier", "Predicted as", "with probability")
	    HTML(res_final, file = "nifpred.html", align = "center", Border = 1, innerBorder = 1)
		}
    	
	}
	}
	###########end#########
	
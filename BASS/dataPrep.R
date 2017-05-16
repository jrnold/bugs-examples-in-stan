## read the High School and beyond data

data1 <- scan(file="HSB1.DAT",
              what=list(school=0,
                minority=0,
                female=0,
                ses=0,
                math=0))
data1 <- as.data.frame(data1)
summary(data1)

data2 <- scan(file="HSB2.DAT",
              what=list(school=0,
                size=0,
                cath=0,
                academic=0,
                disclimate=0,
                minority=0,
                meanses=0))
data2 <- as.data.frame(data2)

using Catalyst
Retina = @reaction_network begin
    µ, 0 --> Sci
    φ, 0 --> Sci1
    β, 0 --> Sci2 
    θ, Ce --> Ce + C 
    λ, Ce1 --> Ce1 + C1
    κ, Ce2 --> Ce2 + C2 
    η, Sci + C --> 0
    ψ, Sci1 + C1 --> 0
    ɛ, Sci2 + C2 --> 0
    k1, Sci --> Sci + Sr + Sh + Sp + Sb 
    p1, Sci1 --> Sci1 + Sb1 + Sh1 
    p26, Sci2 --> Sci2 + Sh2 + Sr2 + Sp2
    k17, Sr --> Sr + R 
    p33, Sr2 --> Sr2 + R2
    k18, Sh --> Sh + HR 
    p7, Sh1 --> Sh1 + HR1
    p19, Sh2 --> Sh2 + HR2
    k19, Sp --> Sp + P 
    p37, Sp2 --> Sp2 + P2
    k28, Sb --> Sb + Rb  
    p9, Sb1 --> Sb1 + Rb1 
    k2, Sr --> 0
    p32, Sr2 --> 0
    k10, Sh --> 0
    p3, Sh1 --> 0
    p18, Sh2 --> 0
    k13, Sp --> 0
    p31, Sp2 --> 0
    k27, Sb --> 0
    p2, Sb1 --> 0
    k3, R + CL --> R + Cf 
    p34, R2 + CL2 --> R2 + Cf2 
    k26, Rb + CH --> Rb + Cf 
    p25, Rb1 + CH1 --> Rb1 + Cf1 
    k11, R --> 0
    p35, R2 --> 0
    k29, Rb --> 0
    p10, Rb1 --> 0
    k14, P + R --> P 
    k16, P --> 0
    p39, P2 + R2 --> P2 
    p38, P2 --> 0
    k4, Cf --> Cp 
    p4, Cf1 --> Cp1 
    p36, Cf2 --> Cp2 
    p11, Cp1 --> CO 
    (k5, k6), Cp <--> Ce 
    (p5, p6), Cp1 <--> Ce1 
    (p28, p29), Cp2 <--> Ce2 
    (k7, k8), Ce <--> E 
    (p16, p17), Ce1 <--> E1 
    α, 0 --> H 
    ω, 0 --> H1 
    δ, 0 --> H2 
    k9, H + HR --> HR + Ce 
    p12, H1 + HR1 --> HR1 + Ce1 
    p21, H2 + HR2 --> HR2 + Ce2 
    k15, Ce + HR --> Ce 
    p20, Ce2 + HR2 --> Ce2
    k12, Cp + CO --> Cg 
    k21, Ce --> 0
    p15, Ce1 --> 0
    p13, Ce2 --> 0
    k22, Ce + Cg --> B 
    k24, B --> 0
end
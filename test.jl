aux_lamb = 2.0 .^ collect(1:21);
aux_dist = 0.1 .^ collect(1:21);

ret1 = 0;
ret1 += XCOrderD1(aux_lamb,aux_dist) / factorial(big(3));
ret1 -= XCOrderD2(aux_lamb,aux_dist) / factorial(big(5));
ret1 += XCOrderD3(aux_lamb,aux_dist) / factorial(big(7));
ret1 -= XCOrderD4(aux_lamb,aux_dist) / factorial(big(9));
ret1 += XCOrderD5(aux_lamb,aux_dist) / factorial(big(11));
ret1 -= XCOrderD6(aux_lamb,aux_dist) / factorial(big(13));
ret1 += XCOrderD7(aux_lamb,aux_dist) / factorial(big(15));
ret1 -= XCOrderD8(aux_lamb,aux_dist) / factorial(big(17));
ret1 += XCOrderD9(aux_lamb,aux_dist) / factorial(big(19));
ret1 -= XCOrderD10(aux_lamb,aux_dist) / factorial(big(21));
ret1 *= exp(-aux_lamb[1] * aux_dist[2]);

ret2 = 0;
ret2 += XCOrder1(aux_lamb) / factorial(big(2));
ret2 -= XCOrder2(aux_lamb,aux_dist) / factorial(big(4));
ret2 += XCOrder3(aux_lamb,aux_dist) / factorial(big(6));
ret2 -= XCOrder4(aux_lamb,aux_dist) / factorial(big(8));
ret2 += XCOrder5(aux_lamb,aux_dist) / factorial(big(10));
ret2 -= XCOrder6(aux_lamb,aux_dist) / factorial(big(12));
ret2 += XCOrder7(aux_lamb,aux_dist) / factorial(big(14));
ret2 -= XCOrder8(aux_lamb,aux_dist) / factorial(big(16));
ret2 += XCOrder9(aux_lamb,aux_dist) / factorial(big(18));
ret2 -= XCOrder10(aux_lamb,aux_dist) / factorial(big(20));
ret2 *= exp(-aux_lamb[1] * aux_dist[2]);


model1_xc_coeffs_ECP = ReadXCCoeffs("ECP_XCCoeffs.txt");
model1_xc_coeffs_FullE = ReadXCCoeffs("FullE_XCCoeffs.txt");
model1_xc_coeffs_ECP_Pol = ReadXCCoeffs("ECP_XCCoeffs_Pol.txt");
model1_xc_coeffs_FullE_Pol = ReadXCCoeffs("FullE_XCCoeffs_Pol.txt");

# model1_xc_coeffs_ECP[1:11] .* gamma.(1.0 .+ 2.0 .* collect(0:10))
# model1_xc_coeffs_ECP[12:22] .* gamma.(1.0 .+ 2.0 .* collect(0:10))
# model1_xc_coeffs_ECP[23:33] .* gamma.(2.0 .+ 2.0 .* collect(0:10))
# model1_xc_coeffs_ECP[34:44] .* gamma.(2.0 .+ 2.0 .* collect(0:10))

# model1_xc_coeffs_ECP_Pol[2] .* 2
# model1_xc_coeffs_ECP_Pol[13] .* 2
# model1_xc_coeffs_ECP_Pol[24] .* 6
# model1_xc_coeffs_ECP_Pol[35] .* 6

# model1_xc_coeffs_ECP[2] .* 2
# model1_xc_coeffs_ECP[13] .* 2
# model1_xc_coeffs_ECP[24] .* 6
# model1_xc_coeffs_ECP[35] .* 6

@printf "%10.6E   " (-model1_xc_coeffs_ECP_Pol[2] .* 2);
@printf "%10.6E   " (-model1_xc_coeffs_ECP_Pol[13] .* 2);
@printf "%10.6E   " (model1_xc_coeffs_ECP_Pol[24] .* 6);
@printf "%10.6E   " (model1_xc_coeffs_ECP_Pol[35] .* 6);

@printf "%10.6E   " (-model1_xc_coeffs_ECP[2] .* 2);
@printf "%10.6E   " (-model1_xc_coeffs_ECP[13] .* 2);
@printf "%10.6E   " (model1_xc_coeffs_ECP[24] .* 6);
@printf "%10.6E \n" (model1_xc_coeffs_ECP[35] .* 6);


@printf "%10.6E   " (-model1_xc_coeffs_FullE_Pol[2] .* 2);
@printf "%10.6E   " (-model1_xc_coeffs_FullE_Pol[13] .* 2);
@printf "%10.6E   " (model1_xc_coeffs_FullE_Pol[24] .* 6);
@printf "%10.6E   " (model1_xc_coeffs_FullE_Pol[35] .* 6);

@printf "%10.6E   " (-model1_xc_coeffs_FullE[2] .* 2);
@printf "%10.6E   " (-model1_xc_coeffs_FullE[13] .* 2);
@printf "%10.6E   " (model1_xc_coeffs_FullE[24] .* 6);
@printf "%10.6E \n" (model1_xc_coeffs_FullE[35] .* 6);

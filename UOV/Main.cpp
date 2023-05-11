#include<NTL/ZZ.h> //pre datovy typ ZZ
#include<NTL/ZZ_p.h>
#include<NTL/vec_ZZ_p.h>
#include<NTL/mat_ZZ_p.h>
#include<NTL/GF2.h>
#include<NTL/vec_GF2.h>
#include<NTL/mat_GF2.h>
#include <iostream>
#include <tuple>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>
#include <NTL/vec_vec_GF2.h>
#include "riesenie_sustavy.h"
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include "UOV_pert.h"
#include "UOV.h"
#include "message_hash.h"
#include "sha256.h"
#include "sha512.h"
#include <vector>

using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib> // for exit function
#include <time.h> //for time measurement
using namespace std;

NTL_CLIENT

int main(int argc, char* argv[]) {
	GF2X modulus;
	long mod = 6;
	long v = 56; long o = 48; long t = 3;
	publicKey_p pk;
	privateKey_p sk;

	publicKey puk;
	privateKey prk;

	Vec<GF2E> dokument;
	Vec<GF2E> podpis;

	if (argc == 1) {
		cout << "No extra Command Line Argument passed "
			"other than program name"
			<< endl;
	}
	if (argc >= 5) {
		std::string mode = argv[1];

		if (mode == "keygen" && argc >= 5) {
			// Parse the parameters
			mod = std::stoi(argv[2]);
			v = std::stoi(argv[3]);
			o = std::stoi(argv[4]);

			// Initialize modulus
			GF2X modulus;
			BuildIrred(modulus, mod);
			GF2E::init(modulus);

			// Execute KeyGen function
			KeyGen(puk, prk, o, v);
			std::ofstream publicKeyFile("key.puk", std::ios::binary);
			if (publicKeyFile.is_open()) {
				publicKeyFile << puk.Q;
				publicKeyFile << puk.L;
				publicKeyFile << puk.A;
				std::cout << "Public key written to key.puk" << std::endl;
			}
			else {
				std::cout << "Unable to open key.puk for writing." << std::endl;
			}
			std::ofstream privateKeyFile("key.priv", std::ios::binary);
			if (privateKeyFile.is_open()) {
				privateKeyFile << prk.Q;
				privateKeyFile << prk.L;
				privateKeyFile << prk.A;
				privateKeyFile << prk.A_T;
				privateKeyFile << prk.b_T;
				privateKeyFile << prk.A_S;
				privateKeyFile << prk.b_S;
				privateKeyFile.close();
				std::cout << "Private key saved to " << "key.priv" << std::endl;
			}
			else {
				std::cout << "Unable to open " << "key.priv" << " for writing." << std::endl;
			}
		}
		else if (mode == "keygenp" && argc >= 6) {
			// Parse the parameters
			mod = std::stoi(argv[2]);
			int v = std::stoi(argv[3]);
			int o = std::stoi(argv[4]);
			int t = std::stoi(argv[5]);

			// Initialize modulus
			GF2X modulus;
			BuildIrred(modulus, mod);
			GF2E::init(modulus);

			// Execute KeyGen_p function
			KeyGen_p(pk, sk, o, v, t);
			std::ofstream publicKeyFile("key.puk", std::ios::binary);
			if (publicKeyFile.is_open()) {
				publicKeyFile << pk.Q;
				publicKeyFile << pk.L;
				publicKeyFile << pk.A;
				std::cout << "Public key written to key.puk" << std::endl;
			}
			else {
				std::cout << "Unable to open key.puk for writing." << std::endl;
			}
			std::ofstream privateKeyFile("key.priv", std::ios::binary);
			if (privateKeyFile.is_open()) {
				privateKeyFile << sk.Q;
				privateKeyFile << sk.Q_wo_z;
				privateKeyFile << sk.polynomy_z;
				privateKeyFile << sk.lambdas;
				privateKeyFile << sk.L;
				privateKeyFile << sk.A;
				privateKeyFile << sk.A_T;
				privateKeyFile << sk.b_T;
				privateKeyFile << sk.A_S;
				privateKeyFile << sk.b_S;
				privateKeyFile.close();
				std::cout << "Private key saved to " << "key.priv" << std::endl;
			}
			else {
				std::cout << "Unable to open " << "key.priv" << " for writing." << std::endl;
			}
		}
		else if (mode == "sign" && argc >= 7) {
			mod = std::stoi(argv[2]);
			v = std::stoi(argv[3]);
			o = std::stoi(argv[4]);

			// Initialize modulus
			GF2X modulus;
			BuildIrred(modulus, mod);
			GF2E::init(modulus);

			hash_file512(dokument, argv[6], o, mod);
			std::ifstream privateKeyFile(argv[5], std::ios::binary);
			if (privateKeyFile.is_open()) {
				privateKeyFile >> prk.Q;
				privateKeyFile >> prk.L;
				privateKeyFile >> prk.A;
				privateKeyFile >> prk.A_T;
				privateKeyFile >> prk.b_T;
				privateKeyFile >> prk.A_S;
				privateKeyFile >> prk.b_S;
				privateKeyFile.close();
			}
			sign(podpis, prk, dokument, v, o);

			std::ofstream signFile("sign.uov", std::ios::binary);
			if (signFile.is_open()) {
				signFile << podpis;
				signFile.close();
				std::cout << "Sign written to sign.uov" << std::endl;
			}
			else {
				std::cout << "Unable to open file for writing." << std::endl;
			}
		}
		else if (mode == "signp" && argc >= 8) {
			mod = std::stoi(argv[2]);
			v = std::stoi(argv[3]);
			o = std::stoi(argv[4]);
			t = std::stoi(argv[5]);

			// Initialize modulus
			GF2X modulus;
			BuildIrred(modulus, mod);
			GF2E::init(modulus);

			hash_file512(dokument, argv[7], o, mod);
			std::ifstream privateKeyFile(argv[6], std::ios::binary);
			if (privateKeyFile.is_open()) {
				privateKeyFile >> sk.Q;
				privateKeyFile >> sk.Q_wo_z;
				privateKeyFile >> sk.polynomy_z;
				privateKeyFile >> sk.lambdas;
				privateKeyFile >> sk.L;
				privateKeyFile >> sk.A;
				privateKeyFile >> sk.A_T;
				privateKeyFile >> sk.b_T;
				privateKeyFile >> sk.A_S;
				privateKeyFile >> sk.b_S;
				privateKeyFile.close();
				std::cout << "Private key loaded from " << argv[6] << std::endl;
			}
			else {
				std::cout << "Unable to open " << argv[6] << " for reading." << std::endl;
			}

			sign_p_v2(podpis, sk, dokument, v, o, t);

			std::ofstream signFile("sign.uovp", std::ios::binary);
			if (signFile.is_open()) {
				signFile << podpis;
				signFile.close();
				std::cout << "Sign written to sign.uovp" << std::endl;
			}
			else {
				std::cout << "Unable to open file for writing." << std::endl;
			}
		}
		else if (mode == "verify" && argc >= 7) {
			mod = std::stoi(argv[2]);
			o = std::stoi(argv[3]);

			// Initialize modulus
			GF2X modulus;
			BuildIrred(modulus, mod);
			GF2E::init(modulus);

			hash_file512(dokument, argv[5], o, mod);

			std::ifstream publicKeyFile(argv[6], std::ios::binary);
			if (publicKeyFile.is_open()) {
				publicKeyFile >> puk.Q;
				publicKeyFile >> puk.L;
				publicKeyFile >> puk.A;
				publicKeyFile.close();

			}
			std::ifstream signFile(argv[4], std::ios::binary);
			if (signFile.is_open()) {
				signFile >> podpis;
				signFile.close();

			}
			verify(podpis, dokument, puk, o);
		}
		else if (mode == "verifyp" && argc >= 7) {
			mod = std::stoi(argv[2]);
			o = std::stoi(argv[3]);

			// Initialize modulus
			GF2X modulus;
			BuildIrred(modulus, mod);
			GF2E::init(modulus);
			hash_file512(dokument, argv[5], o, mod);

			std::ifstream publicKeyFile(argv[6], std::ios::binary);
			if (publicKeyFile.is_open()) {
				publicKeyFile >> pk.Q;
				publicKeyFile >> pk.L;
				publicKeyFile >> pk.A;
				publicKeyFile.close();

			}
			std::ifstream signFile(argv[4], std::ios::binary);
			if (signFile.is_open()) {
				signFile >> podpis;
				signFile.close();

			}
			verify_p(podpis, dokument, pk, o);
		}
		else if (mode == "all" && argc >= 7) {
			mod = std::stoi(argv[2]);
			int v = std::stoi(argv[3]);
			int o = std::stoi(argv[4]);
			int t = std::stoi(argv[5]);

			// Initialize modulus
			GF2X modulus;
			BuildIrred(modulus, mod);
			GF2E::init(modulus);

			// Execute KeyGen_p function
			std::cout << "Generating keys ..." << std::endl;
			KeyGen_p(pk, sk, o, v, t);
			std::cout << "hashing file and genereting sign" << std::endl;
			hash_file512(dokument, argv[6], o, mod);
			sign_p_v2(podpis, sk, dokument, v, o, t);
			std::cout << "Sign done" << std::endl;
			verify_p(podpis, dokument, pk, o);

		}

		else {
			std::cout << "Invalid mode or insufficient parameters or other error" << std::endl;
		}

	}
}

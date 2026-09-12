// Author: John Cesar
// Date: 16/09/2022
// Description: Program to read in a user-inputed TPPT coincidence file and
// split the data by ChipID and also skip lines with energies below a chosen
// threshold of 10 DAQ units

// HOW TO USE: open a terminal in the directory which contains both the _coinc.dat 
// file to be pre-processed as well as this macro. Make sure that the variable 
// below "std::string filename" is initialized to the name of the file to be 
// pre-processed. Save the macro and then type, in the terminal:
//   root -l Data_Preprocess_v2.C

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

int Data_Preprocess_v2(){

  printf("...Initializations...\n");
  // Initialize variables for file reading
  int Size1, Size2, Index1, Index2,
    AbsID1, AbsID2, SlaveID1, SlaveID2, ChipID1, ChipID2, ChanID1, ChanID2,
    PCBChanID1, PCBChanID2, GeoID1, GeoID2;
  long long Time1, Time2;
  double Energy1, Energy2, E_Threshold;

  // This is the variable that defines the basic energy pre-cut
  // We use to to skip any lines where either E1 or E2 are less than E_threshold
  E_Threshold = 10.0;

  // Initialize variables for input and output streams
  string line;
  std::ifstream in;
  std::ofstream out1, out2, out3, out4, out5, out6, out7, out8, out9, out10,
    out11, out12, out13, out14, out15, out16, out17, out18, out19, out20,
    out21, out22, out23, out24, out25, out26, out27, out28, out29, out30,
    out31, out32, out33, out34, out35, out36, out37, out38, out39, out40,
    out41, out42, out43, out44, out45, out46, out47, out48;

  // Initialize some variable useful for debugging
  bool Debug = 0;
  int lines_skipped = 0, lines_total = 0, lines_cut = 0;
  //  std::string filename = "Test_coinc.dat";
  std::string filename = "HWTriggerTest8.2FullScanner10min45_coinc.dat";

  // Prompt the user for the input filename and read in to a string
  // std::cout << "While file would you like to pre-process?\n";
  // std::cout << "Input full filename: ";
  // std::cin >> filename;

  // Open the input file stream
  in.open(filename);

  // Perform some string manipulation to truncate the "_coinc.dat"
  std::size_t split = filename.find_last_of("."); // find the last "period" in the filename
  std::string truncated_filename = filename.substr(0,split-6); // new string goes from 0th Index
                                                               // to 6 spaces before the last period
                                                               // in order to truncate "_coinc"
  if (Debug) printf("...Data will output to filename starting with: %s...\n",truncated_filename.c_str());

  // Open the output file streams - Slave0 outputs
  out1.open(truncated_filename+"_Slave0Chip0.dat");
  out2.open(truncated_filename+"_Slave0Chip1.dat");
  out3.open(truncated_filename+"_Slave0Chip2.dat");
  out4.open(truncated_filename+"_Slave0Chip3.dat");
  out5.open(truncated_filename+"_Slave0Chip4.dat");
  out6.open(truncated_filename+"_Slave0Chip5.dat");
  out7.open(truncated_filename+"_Slave0Chip6.dat");
  out8.open(truncated_filename+"_Slave0Chip7.dat");
  out9.open(truncated_filename+"_Slave0Chip8.dat");
  out10.open(truncated_filename+"_Slave0Chip9.dat");
  out11.open(truncated_filename+"_Slave0Chip10.dat");
  out12.open(truncated_filename+"_Slave0Chip11.dat");
  out13.open(truncated_filename+"_Slave0Chip12.dat");
  out14.open(truncated_filename+"_Slave0Chip13.dat");
  out15.open(truncated_filename+"_Slave0Chip14.dat");
  out16.open(truncated_filename+"_Slave0Chip15.dat");

  // Open the output file streams - Slave1 outputs
  out17.open(truncated_filename+"_Slave1Chip0.dat");
  out18.open(truncated_filename+"_Slave1Chip1.dat");
  out19.open(truncated_filename+"_Slave1Chip2.dat");
  out20.open(truncated_filename+"_Slave1Chip3.dat");
  out21.open(truncated_filename+"_Slave1Chip4.dat");
  out22.open(truncated_filename+"_Slave1Chip5.dat");
  out23.open(truncated_filename+"_Slave1Chip6.dat");
  out24.open(truncated_filename+"_Slave1Chip7.dat");
  out25.open(truncated_filename+"_Slave1Chip8.dat");
  out26.open(truncated_filename+"_Slave1Chip9.dat");
  out27.open(truncated_filename+"_Slave1Chip10.dat");
  out28.open(truncated_filename+"_Slave1Chip11.dat");
  out29.open(truncated_filename+"_Slave1Chip12.dat");
  out30.open(truncated_filename+"_Slave1Chip13.dat");
  out31.open(truncated_filename+"_Slave1Chip14.dat");
  out32.open(truncated_filename+"_Slave1Chip15.dat");

  // Open the output file streams - Slave2 outputs
  out33.open(truncated_filename+"_Slave2Chip0.dat");
  out34.open(truncated_filename+"_Slave2Chip1.dat");
  out35.open(truncated_filename+"_Slave2Chip2.dat");
  out36.open(truncated_filename+"_Slave2Chip3.dat");
  out37.open(truncated_filename+"_Slave2Chip4.dat");
  out38.open(truncated_filename+"_Slave2Chip5.dat");
  out39.open(truncated_filename+"_Slave2Chip6.dat");
  out40.open(truncated_filename+"_Slave2Chip7.dat");
  out41.open(truncated_filename+"_Slave2Chip8.dat");
  out42.open(truncated_filename+"_Slave2Chip9.dat");
  out43.open(truncated_filename+"_Slave2Chip10.dat");
  out44.open(truncated_filename+"_Slave2Chip11.dat");
  out45.open(truncated_filename+"_Slave2Chip12.dat");
  out46.open(truncated_filename+"_Slave2Chip13.dat");
  out47.open(truncated_filename+"_Slave2Chip14.dat");
  out48.open(truncated_filename+"_Slave2Chip15.dat");

  printf("...Reading in data file...\n");
  // Read in *_coinc.dat file, the format for the input file is:
  // col0        col1         col2    col3      col4     col5       col6         col7    col8      col9
  // Size_Grp1   Index_Grp1   Time1   Energy1   AbsID1   Size_Grp2  Index_Grp2   Time2   Energy2   AbsID2
  // while(in >> Size1 >> Index1 >> Time1 >> Energy1 >> AbsID1 >> Size2 >> Index2 >> Time2 >> Energy2 >> AbsID2){
  while(getline(in,line)){

    // Count the number of total lines to compare to the skipped lines
    lines_total++;

    // Skip line if there is an error reading in one of the columns
    istringstream iss(line);
    if (!(iss >> Size1 >> Index1 >> Time1 >> Energy1 >> AbsID1 >> Size2 >> Index2 >> Time2 >> Energy2 >> AbsID2)){
            lines_skipped++;
            continue;	
    }	

    // Cut line if either energy is below the threshold defined above
    if((Energy1 < E_Threshold) || (Energy2 < E_Threshold)){
      if(Debug) printf("  ...Skipping line with E1 or E2 < E_threshold: E1 = %f | E2 = %f\n",Energy1,Energy2);
      // Count the number of skipped lines
      lines_cut++;
      continue;
    }
    // if the energy cut is passed, use the AbsID2 (from the "right-side crescent")
    // to determine which Slave and Chip the event belongs to and send to the appropriate
    // output file
    else{
      if(Debug) printf("  ...Energy cut passed...\n");

      // Slave0 outputs - AbsID within 0-1024
      if( (AbsID2 >= 0) & (AbsID2 < 64) ) out1 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 64) & (AbsID2 < 128) ) out2 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 128) & (AbsID2 < 192) ) out3 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 192) & (AbsID2 < 256) ) out4 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 256) & (AbsID2 < 320) ) out5 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 320) & (AbsID2 < 384) ) out6 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 384) & (AbsID2 < 448) ) out7 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 448) & (AbsID2 < 512) ) out8 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 512) & (AbsID2 < 576) ) out9 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 576) & (AbsID2 < 640) ) out10 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 640) & (AbsID2 < 704) ) out11 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 704) & (AbsID2 < 768) ) out12 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 768) & (AbsID2 < 832) ) out13 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 832) & (AbsID2 < 896) ) out14 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 896) & (AbsID2 < 960) ) out15 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 960) & (AbsID2 < 1024) ) out16 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";

      // Slave1 outputs - AbsID within 4096-5120
      else if( (AbsID2 >= 4096) & (AbsID2 < 4160) ) out17 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4160) & (AbsID2 < 4224) ) out18 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4224) & (AbsID2 < 4288) ) out19 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4288) & (AbsID2 < 4352) ) out20 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4352) & (AbsID2 < 4416) ) out21 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4416) & (AbsID2 < 4480) ) out22 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4480) & (AbsID2 < 4544) ) out23 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4544) & (AbsID2 < 4608) ) out24 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4608) & (AbsID2 < 4672) ) out25 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4672) & (AbsID2 < 4736) ) out26 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4736) & (AbsID2 < 4800) ) out27 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4800) & (AbsID2 < 4864) ) out28 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4864) & (AbsID2 < 4928) ) out29 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4928) & (AbsID2 < 4992) ) out30 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 4992) & (AbsID2 < 5056) ) out31 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 5056) & (AbsID2 < 5120) ) out32 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";

      // Slave2 outputs - AbsID within 8192-9216
      else if( (AbsID2 >= 8192) & (AbsID2 < 8256) ) out33 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8256) & (AbsID2 < 8320) ) out34 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8320) & (AbsID2 < 8384) ) out35 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8384) & (AbsID2 < 8448) ) out36 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8448) & (AbsID2 < 8512) ) out37 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8512) & (AbsID2 < 8576) ) out38 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8576) & (AbsID2 < 8640) ) out39 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8640) & (AbsID2 < 8704) ) out40 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8704) & (AbsID2 < 8768) ) out41 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8768) & (AbsID2 < 8832) ) out42 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8832) & (AbsID2 < 8896) ) out43 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8896) & (AbsID2 < 8960) ) out44 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 8960) & (AbsID2 < 9024) ) out45 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 9024) & (AbsID2 < 9088) ) out46 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 9088) & (AbsID2 < 9152) ) out47 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";
      else if( (AbsID2 >= 9152) & (AbsID2 < 9216) ) out48 << Size1 << "\t" << Index1 << "\t" << Time1 << "\t" << Energy1 << "\t" << AbsID1 << "\t" << Size2 << "\t" << Index2 << "\t" << Time2 << "\t" << Energy2 << "\t" << AbsID2 << "\n";

      else continue;
    }
  }
  // Close all the input and output streams
  in.close();
  out1.close();
  out2.close();
  out3.close();
  out4.close();
  out5.close();
  out6.close();
  out7.close();
  out8.close();
  out9.close();
  out10.close();
  out11.close();
  out12.close();
  out13.close();
  out14.close();
  out15.close();
  out16.close();
  out17.close();
  out18.close();
  out19.close();
  out20.close();
  out21.close();
  out22.close();
  out23.close();
  out24.close();
  out25.close();
  out26.close();
  out27.close();
  out28.close();
  out29.close();
  out30.close();
  out31.close();
  out32.close();
  out33.close();
  out34.close();
  out35.close();
  out36.close();
  out37.close();
  out38.close();
  out39.close();
  out40.close();
  out41.close();
  out42.close();
  out43.close();
  out44.close();
  out45.close();
  out46.close();
  out47.close();
  out48.close();

  // Print a summary of the total number of lines that were skipped fro mthe E_threshold cut
  printf("\n\n There were %d total lines...",lines_total);
  printf("\n   ...of which %d were skipped due to errors...",lines_skipped);
  printf("\n   ...and of which %d were cut due to low energies...",lines_cut);
  printf("\n\n A total of %f percent of lines were skipped and %f percent were cut \n\n",double(lines_skipped)/double(lines_total)*100,double(lines_cut)/double(lines_total)*100);

  return(23);
}

#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <zlib.h>
// gzfile
// Thanks https://stackoverflow.com/questions/624250/how-do-i-read-write-gzipped-files-in-c
// And this solution https://gist.github.com/piti118/1508048
#include "gzstream.h"

using namespace Rcpp;

// cf stackoverflow.com/questions/874134
bool hasEnding(std::string const &full, std::string const &ending) {
  // Determine if we have a .gz file ending
    if (full.length() >= ending.length()) {
        return(0 == full.compare(full.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

// ------------------------------------------------------------------
// [[Rcpp::export]]
List readfilecpp(std::string path, int nhaps, int nsnps,double scale,bool verbose)
{
  // Data to be computed
  // SNP information
  NumericVector snpsums (nsnps);
  NumericVector snpsumsq (nsnps);
  StringVector snpnames (nsnps);
  // Haplotype information
  NumericVector hapsums (nhaps);
  NumericVector hapsumsq (nhaps);
  StringVector hapnames (nhaps);

  // Two file streams. We make a pointer to the right one
  std::ifstream myfileraw(path.c_str());
  igzstream in(path.c_str());

  // Decide whether to treat the file as gzip or not
  std::istream *myfile;
  if (hasEnding(path, ".gz")) {
    if(verbose) Rprintf("File %s is gzipped...\n",path.c_str());
    myfile = &in;//new std::istream(&in); // access the file via the un-gzip stream
  }else{
    if(verbose) Rprintf("File %s is plain text...\n",path.c_str());
    myfile = &myfileraw; // access the file directly
  }
  // Now myfile can be used invisbly to which stream is used.
  
  if(verbose) Rprintf("Reading header...\n");
  // Read the header
  std::string line;
  std::getline(*myfile,line,'\n'); // Get the header
  std::stringstream iss0(line); // take the line into a stringstream
  std::string val;
  std::getline(iss0,val,' '); //Omit the "top left" value, name for the row names column
  for (int col=0; col<nsnps; ++col )
    {
      std::string val;
      std::getline(iss0,val,' '); //reads the stringstream line and separate by spaces
      snpnames[nsnps - col - 1]=val;
    }
  // Read the data
  for (int row=0; row<nhaps; ++row)
    {
      if(verbose) Rprintf("row %d of %d...\n",row,nhaps);
      std::string line;
      std::getline(*myfile,line,'\n'); //Starts at the second line
      
      if(!myfile->good() ) {//If end of rows then break
	Rprintf("Error with file %s...\n",path.c_str());
	break;
      }
      std::stringstream iss(line); // take the line into a stringstream
      std::string val;
      std::getline(iss,val,' '); ///skips the first column (row names)
      hapnames[row]=val;
      
      for (int col=0; col<nsnps; ++col )
	{
	  std::string val;
	  std::getline(iss,val,' '); //reads the stringstream line and separate by spaces

	  std::stringstream convertor(val); //get the results into another stringstream convertor
	  double dval;
	  convertor >>dval;
	  hapsums[row] += dval * scale;
	  hapsumsq[row] += (dval * scale)*(dval * scale);
	  snpsums[nsnps - col - 1] += dval * scale;
	  snpsumsq[nsnps - col - 1] += (dval * scale)*(dval * scale);
	}
    }
  return( List::create(Named("snpsums") = snpsums,
		       Named("hapsums") = hapsums,
		       Named("snpsumsq") = snpsumsq,
		       Named("hapsumsq") = hapsumsq,
		       Named("snpnames") = snpnames,
		       Named("hapnames") = hapnames
		       ));
}

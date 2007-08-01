/*****************************************************************************/
void printToFile(TH1F *h1, const char* fileName)
{
  ofstream file(fileName);
  for(int i = 1; i <= h1->GetNbinsX(); i++)
    file << " " << h1->GetBinCenter(i)
         << " " << h1->GetBinContent(i)
         << " " << h1->GetBinError(i)
         << endl;
  file.close();
}

/*****************************************************************************/
void printToFile(TH2F *h2, const char* fileName)
{
  ofstream file(fileName);
  for(int i = 1; i <= h2->GetNbinsX(); i++)
  {
    for(int j = 1; j <= h2->GetNbinsY(); j++)
      file << " " << h2->GetXaxis()->GetBinCenter(i)
           << " " << h2->GetYaxis()->GetBinCenter(j)
           << " " << h2->GetBinContent(i,j)
           << " " << h2->GetBinError(i,j)
           << endl;
    file << endl;
  }
  file.close();
}

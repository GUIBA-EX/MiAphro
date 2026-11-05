// Independent exhaustive all-window scanner: no tail prefilter, no Python scanner imports.
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
using namespace std;
int mask(char c){switch(c){case'A':return 1;case'C':return 2;case'G':return 4;case'T':return 8;case'R':return 5;case'Y':return 10;case'S':return 6;case'W':return 9;case'K':return 12;case'M':return 3;case'B':return 14;case'D':return 13;case'H':return 11;case'V':return 7;default:return 15;}}
int comp(int x){return ((x&1)?8:0)|((x&2)?4:0)|((x&4)?2:0)|((x&8)?1:0);}
struct Rec{string id;vector<int>s;};
vector<Rec> fasta(string f){ifstream in(f);string l;vector<Rec> r;while(getline(in,l)){if(l.empty())continue;if(l[0]=='>'){istringstream ss(l.substr(1));string id;ss>>id;r.push_back({id,{}});}else for(char c:l)if(!isspace(c))r.back().s.push_back(mask(toupper(c)));}return r;}
struct Hit{int pos,mm,last3,amb;};
vector<int> primer(string q,bool reverse,bool gt){vector<int> m;for(char c:q)m.push_back(mask(c));if(reverse){std::reverse(m.begin(),m.end());for(auto &x:m)x=comp(x);}if(gt)for(auto&x:m){int old=x;if(!reverse){if(old&4)x|=1;if(old&8)x|=2;}else{if(old&2)x|=8;if(old&1)x|=4;}}return m;}
vector<Hit> scan(const vector<int>&s,const vector<int>&q,bool reverse,int limit){vector<Hit> h;for(int i=0;i+(int)q.size()<=(int)s.size();i++){int mm=0,m3=0,amb=0;for(int j=0;j<(int)q.size();j++){int x=s[i+j];amb+=x!=1&&x!=2&&x!=4&&x!=8;int mis=(x&q[j])==0;mm+=mis;if(reverse?j<3:j>=(int)q.size()-3)m3+=mis;if(mm>limit)break;}if(mm<=limit)h.push_back({i,mm,m3,amb});}return h;}
struct P{string name,F,R;};
int main(int argc,char**argv){string root=argv[1];ifstream in(root+"/independent_primers.tsv");string l;vector<P> ps;while(getline(in,l)){istringstream ss(l);P p;ss>>p.name>>p.F>>p.R;ps.push_back(p);}ofstream barcode(root+"/target_barcodes.tsv");barcode<<"name\taccession\tsequence\n";ofstream out(root+"/independent_hits.tsv");out<<"name\trole\tmodel\taccession\texact\tle2last3\tle2\tle3\tle4\tambiguous\tsizes\n";
for(string role:{"target","background"}){auto records=fasta(root+"/"+role+".fasta");for(auto&p:ps)for(int gt=0;gt<(role=="target"?1:2);gt++){auto F=primer(p.F,false,gt),R=primer(p.R,true,gt);int limit=role=="target"||gt?2:4;for(auto &rec:records){bool exact=false,m3=false,le2=false,le3=false,le4=false,amb=false;set<int> sizes;set<string> bars;for(int ori=0;ori<2;ori++){vector<int>s=rec.s;if(ori){reverse(s.begin(),s.end());for(auto&x:s)x=comp(x);}auto fh=scan(s,F,false,limit),rh=scan(s,R,true,limit);for(auto&f:fh)for(auto&r:rh){int inner=r.pos-f.pos-(int)F.size(),total=r.pos+(int)R.size()-f.pos;if(role=="target"?(inner<50||total<100||total>500):(inner<0||inner>1000))continue;if(f.amb||r.amb){amb=true;continue;}exact|=f.mm==0&&r.mm==0;le2|=f.mm<=2&&r.mm<=2;le3|=f.mm<=3&&r.mm<=3;le4=true;if(f.mm<=2&&r.mm<=2&&!f.last3&&!r.last3){m3=true;sizes.insert(total);if(role=="target"){string bar;bool clean=true;for(int bp=f.pos+F.size();bp<r.pos;bp++){int b=s[bp];if(b!=1&&b!=2&&b!=4&&b!=8){clean=false;break;}bar+=(b==1?'A':b==2?'C':b==4?'G':'T');}if(clean)bars.insert(bar);}}}}for(auto &bar:bars)barcode<<p.name<<'\t'<<rec.id<<'\t'<<bar<<'\n';if(exact||le2||le3||le4||amb){out<<p.name<<'\t'<<role<<'\t'<<(gt?"GT":"WC")<<'\t'<<rec.id<<'\t'<<exact<<'\t'<<m3<<'\t'<<le2<<'\t'<<le3<<'\t'<<le4<<'\t'<<amb<<'\t';bool first=true;for(int n:sizes){if(!first)out<<',';out<<n;first=false;}out<<'\n';}}cerr<<role<<' '<<p.name<<' '<<gt<<" done\n";}}
}

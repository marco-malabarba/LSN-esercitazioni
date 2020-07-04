#include "popolazione.h"

using namespace std;

popolazione::popolazione (double *x, double *y, int n_cities, int n_popolazione) {
	if (n_cities<0||n_popolazione<0) {
		cerr<<"Errore"<<endl;
		exit(1);
	}
	else {
	
		//setrandom(rnd,k); //Lo faccio dopo
		m_n_cities=n_cities;
		m_n_individui=n_popolazione;
		m_x_city= new double [m_n_cities];
		m_y_city=new double [m_n_cities];
		m_lunghezze= new double [m_n_individui];
		m_prob=new double [m_n_individui];
		m_individui= new int *[m_n_individui];
		for (int i=0; i<m_n_individui; i++) {
			m_individui[i]=new int [m_n_cities];
		}
		for (int i=0; i<m_n_cities; i++) {
			m_x_city[i]=x[i];
			m_y_city[i]=y[i];
		}

		for (int i=0; i<m_n_individui; i++) {
			m_prob[i]=0.;
			m_lunghezze[i]=0.;
			for (int j=0; j<m_n_cities; j++) {
				m_individui[i][j]=0;
			}
		}
			
	}	
}

popolazione::~popolazione() {
	delete []m_x_city;
	delete []m_y_city;
	delete []m_lunghezze;
	delete []m_prob;
	for (int i=0; i<m_n_individui; i++) {
		delete [] m_individui[i];
	}
	delete [] m_individui;
}

void popolazione::setrndm(int k) {
	setrandom(rnd,k);
}

void popolazione::creaindividui(void) {
	for (int i=0; i<m_n_individui; i++) {
		if (i==0) {
			for (int j=0; j<m_n_cities; j++) {
				m_individui[i][j]=j; //!!!!!!!!!!!
			}
		}
		else {
			for (int j=0; j<m_n_cities; j++) {
				m_individui[i][j]=m_individui[i-1][j];
			}
			int swap_1=(int) rnd.Rannyu(1.,m_n_cities);
			int swap_2=(int) rnd.Rannyu(1.,m_n_cities);
			int temp=m_individui[i][swap_2];
			m_individui[i][swap_2]=m_individui[i][swap_1];
			m_individui[i][swap_1]=temp;
		}
	}
}

void popolazione::stampaindividui() {
	for (int i=0; i<m_n_individui; i++) {
		cout<<i+1<<":"<<'\t';
		for (int j=0; j<m_n_cities; j++) {
			cout<<m_individui[i][j]<<'\t';
		}
		cout<<endl;
	}
	cout<<endl;
}

void popolazione::calcolalunghezze() {
	for (int i=0; i<m_n_individui; i++) {
		m_lunghezze[i]=0.;
		for (int j=0; j<m_n_cities; j++) {
			int p_old, p_new;
			if ((j+1)!=m_n_cities) {
				p_old=m_individui[i][j];
				p_new=m_individui[i][j+1];
				m_lunghezze[i]+=distanza(m_x_city[p_old],m_y_city[p_old],m_x_city[p_new],m_y_city[p_new]);
			}
			else {
				p_old=m_individui[i][j];
				p_new=m_individui[i][0];
				m_lunghezze[i]+=distanza(m_x_city[p_old],m_y_city[p_old],m_x_city[p_new],m_y_city[p_new]);
			}
		}
	}
}	

void popolazione::calcolaprob() {
	int index_min;
	double sum=0;
	for (int i=0; i<m_n_individui; i++) {
		sum+=1./pow(m_lunghezze[i],2);
	}
	for (int i=0; i<m_n_individui; i++) {
		m_prob[i]=1./(sum*pow(m_lunghezze[i],2));
	
	}
}


void popolazione::stampalunghezze() {
	for (int i=0; i<m_n_individui; i++) {
		cout<<i<<" "<<m_lunghezze[i]<<endl;
	}
}

void popolazione::stampaprob() {
	for (int i=0; i<m_n_individui; i++) {
		cout<<i<<" "<<m_prob[i]<<endl;
	}
}

void popolazione::mutazioni() {
	for (int i=0; i<m_n_individui; i++) {
		if (rnd.Rannyu()<p_mut1) {
			int ind_1=(int)rnd.Rannyu(1.,m_n_cities);
			int ind_2=(int)rnd.Rannyu(1.,m_n_cities);
			swap(m_individui[i][ind_1],m_individui[i][ind_2]);
		}
	}
	
	for (int i=0; i<m_n_individui; i++) {
		if (rnd.Rannyu()<p_mut2) {
			int trasl=(int) rnd.Rannyu(1.,m_n_cities);
			int *val=new int [m_n_cities];
			for (int j=0; j<m_n_cities; j++) val[j]=m_individui[i][j];
			for (int j=1; j<m_n_cities; j++) {
				int ind=pbc(j+trasl,m_n_cities);
				m_individui[i][ind]=val[j];
			}
		}
	}

	for (int i=0; i<m_n_individui; i++) {
		if (rnd.Rannyu()<p_mut3) {
			int blocco=(int) rnd.Rannyu(2.,m_n_cities/2.);
			int ind_1=(int) rnd.Rannyu(1.,m_n_cities/2.);
			int ind_2=(int) rnd.Rannyu(m_n_cities/2.,m_n_cities);
			while (ind_2+ind_1>=m_n_cities||ind_1+blocco>=ind_2||ind_2+blocco>=m_n_cities) {
				blocco=(int) rnd.Rannyu(2.,m_n_cities/2.);
				ind_1=(int) rnd.Rannyu(1.,m_n_cities/2.);
				ind_2=(int) rnd.Rannyu(m_n_cities/2.,m_n_cities);
			}
			//cout<<blocco<<" "<<ind_1<<" "<<ind_2<<endl;
			for (int j=0; j<blocco; j++) {
				int temp=m_individui[i][pbc(j+ind_2,m_n_cities)];
				m_individui[i][pbc(j+ind_2,m_n_cities)]=m_individui[i][pbc(j+ind_1,m_n_cities)];
		 		m_individui[i][pbc(j+ind_1,m_n_cities)]=temp;
			}
		}
	}

	for (int i=0; i<m_n_individui; i++) {
		if (rnd.Rannyu()<p_mut4) {
			int n_inv=(int) rnd.Rannyu(2.,m_n_cities);
			int ind=(int) rnd.Rannyu(1.,m_n_cities);
			for (int j=0; j<=(int)(n_inv/2.); j++) {
				swap(m_individui[i][pbc(j+ind,m_n_cities)],m_individui[i][pbc(ind+n_inv-j,m_n_cities)]);
			}
		}
	}
}
	
void popolazione::riproduzione() {
	int **figli=new int *[m_n_individui];
	for (int i=0; i<m_n_individui; i++) {
		figli[i]=new int [m_n_cities];
	}

	for (int i=0; i<m_n_individui; i+=2) {
		double alfa1=rnd.Rannyu();
		double alfa2=rnd.Rannyu();
		int index=0;
		double sum=m_prob[index];
		while (sum<alfa1) {
			index++;
			sum+=m_prob[index];
		}
		int genitore_1=index;

		index=0;
		sum=m_prob[index];
		while (sum<alfa2) {
			index++;
			sum+=m_prob[index];
		}
		int genitore_2=index;
		//cout<<genitore_1+1<<" "<<genitore_2+1<<endl;

		int mezzo=(int)m_n_cities/2;

		for (int j=0; j<mezzo; j++) {
			figli[i][j]=m_individui[genitore_1][j];
			figli[i+1][j]=m_individui[genitore_2][j];
		}
		int p=0;
		for (int j=mezzo; j<m_n_cities; j++) {
			int ind=0;
			while(cerca(m_individui[genitore_2][ind],figli[i],mezzo+p)==true) ind++;
			figli[i][j]=m_individui[genitore_2][ind];
			ind=0;
			while(cerca(m_individui[genitore_1][ind],figli[i+1],mezzo+p)==true) ind++;
			figli[i+1][j]=m_individui[genitore_1][ind];
			p++;			
		}
	}
	for (int i=0; i<m_n_individui; i++) {
		for (int j=0; j<m_n_cities; j++) {
			m_individui[i][j]=figli[i][j];
		}
	}
}

double popolazione::ave_best() {
	double *copy=new double [m_n_individui]; //creo una copia per evitare problemi nelle riproduzioni
	for (int i=0; i<m_n_individui; i++) copy[i]=m_lunghezze[i];
	sort(copy,m_n_individui);
	int mezzo= (int) m_n_individui/2.;
	double ave=0.;
	for (int i=0; i<mezzo; i++) ave+=copy[i];
	ave=ave/(double) mezzo;
	return ave;
}

int popolazione::cerca_best() {
	int index=0;
	for (int i=1; i<m_n_individui; i++) {
		if (m_lunghezze[i]<m_lunghezze[index]) index=i;
	}
	return index;
}

double popolazione::best(){
	int u=cerca_best();
	//cout<<u<<endl;
	return m_lunghezze[u];
	}

void popolazione::stampapercorso (int q, string k) {
	string titolo=(string)"percorso_"+k+(string)".dat";
	ofstream out(titolo);
	for (int i=0; i<m_n_cities; i++) {
		int ind=m_individui[q][i];
		out<<m_x_city[ind]<<" "<<m_y_city[ind]<<endl;
	}
	out<<m_x_city[0]<<" "<<m_y_city[0]<<endl;
	out.close();
}
			
void popolazione::sim_ann(double temp) {
	if (m_n_individui!=1) {
		cerr<<"Errore! Il simulated annealing funziona solo con un singolo individuo!"<<endl;
		exit(1);
	}
	double beta=1./temp;
	calcolalunghezze();
	double l_old=m_lunghezze[0];
	int *perc_old=new int [m_n_cities];
	for (int i=0; i<m_n_cities; i++) perc_old[i]=m_individui[0][i];
	mutazioni();
	calcolalunghezze();
	double l_new=m_lunghezze[0];

	//stepmetropolis(beta,l_new,l_old,m_individui[0],perc_old,m_n_cities,rnd);
	double alfa=fmin(1.,exp(-beta*(pow(l_new,1)-pow(l_old,1))));
	if (rnd.Rannyu()<alfa) {
		//cout<<"*"<<endl;
	}
	else {
		for (int i=0; i<m_n_cities; i++) m_individui[0][i]=perc_old[i];
		calcolalunghezze();
	
	}
	
	
}

void popolazione::set_p_mut (double p1, double p2, double p3, double p4) {
	p_mut1=p1;
	p_mut2=p2;
	p_mut3=p3;
	p_mut4=p4;
}
	
int popolazione::getcities() {
	return m_n_cities;
}

int popolazione::getindividui() {
	return m_n_individui;
}	

int popolazione::getpercorso(int x, int y) {
	return m_individui[x][y];
}

void popolazione::setpercorso(int x, int y, int z) {
	m_individui[x][y]=z;
}

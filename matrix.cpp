#include <iostream>
using namespace std;

class Matrix
{
protected:
    int n;
    int m;
    float* data;
private:
    bool okRow(float* data, int n, int index)
    {
        if(data[index*n+index] != 0)
            return true;
        else
        {
            for(int i=index+1; i<n; i++)
            {
                if(data[i*n+index]!=0)
                {
                    float tmp;
                    for(int u=0; u<n; u++)
                    {
                        tmp = data[index*n+u];
                        if((i - index) % 2 == 1)
                            tmp = -tmp;
                        data[index*n+u] = data[i*n+u];
                        data[i*n+u] = tmp;
                    }
                    return true;
                }
            }
        }
        return false;
    }
    float determinant(float data[], int n)
    {
        for(int i=0; i<n; i++)
        {
            if(okRow(data,n,i))
            {
                for(int j=i+1; j<n; j++)
                {
                    float K = data[j*n+i]/data[i*n+i];
                    for(int k = 0; k<n; k++)
                        data[j*n+k] = data[j*n+k] - K*data[i*n+k];
                }
            }
        }
        float det = 1;
        for(int i=0; i<n; i++)
            det = det*data[i*n+i];
        if(det < 0.0001 && det > -0.0001)
            det = 0;
        return det;
    }

public:
    Matrix()
    {
        this->n = 0;
        this->m = 0;
        this->data = NULL;
    }
    Matrix(int m, int n)
    {
        if( m <=0 || n<=0)
        {
            Matrix();
            return;
        }
        this->m = m;
        this->n = n;
        this->data = new float [n*m];
    }
    Matrix(const Matrix& t) // copy constructor
    {
        this->m = t.m;
        this->n = t.n;
        this->data = new float [m*n];
        for(int i=0; i<m; i++)
            for(int j=0; j<n; j++)
                data[i*n + j] = t.data[i*n + j];
    }
    virtual ~Matrix()
    {
        delete data;
    }
    virtual Matrix operator = (Matrix& t)
    {
        this->m = t.getM();
        this->n = t.getN();
        if(data != NULL)
            delete data;
        this->data = new float [m*n];
        for(int i=0; i<m; i++)
            for(int j=0; j<n; j++)
                data[i*n+j] = t.my_get(i, j);
        return *this;
    }
    virtual Matrix operator+(Matrix& t)
    {
        if(n!=t.getN() || m!=t.getM())
        {
            Matrix sum;
            return sum;
        }
        Matrix sum(m,n);
        for(int i=0; i<n; i++)
            for(int j=0; j<m; j++)
            {
                float z = this->my_get(i,j) + t.my_get(i,j);
                sum.set(i,j,z);
            }
        return sum;
    }
    virtual Matrix operator*(Matrix& t)
    {
        if(this->n != t.getM() || failed() || t.failed())
        {
            Matrix mult;
            return mult;
        }
        Matrix mult(this->m, t.getN());
        for(int m = 0; m < mult.getM(); m++)
            for(int k = 0; k < mult.getN(); k++)
            {
                float sumMult = 0;
                for(int j=0; j<n; j++)
                    sumMult += (my_get(m,j) * t.my_get(j,k));
                mult.set(m,k,sumMult);
            }
        return mult;
    }
    virtual Matrix operator*(float& f)
    {
        if( failed() )
        {
            Matrix mult;
            return mult;
        }
        Matrix mult(m,n);
        for(int i=0; i<m*n; i++)
            mult.data[i] = data[i] * f;
        return mult;
    }
    virtual Matrix operator-(Matrix& t)
    {
        if(n!=t.getN() || m!=t.getM())
        {
            Matrix sum;
            return sum;
        }
        Matrix sum(m,n);
        for(int i=0; i<n; i++)
            for(int j=0; j<m; j++)
            {
                float z = this->my_get(i,j) - t.my_get(i,j);
                sum.set(i,j,z);
            }
        return sum;
    }
    Matrix minor(int I, int J)
    {
        Matrix Minor((n-1),(n-1));
        int count  = 0;
        for(int i=0; i<n*n; i++)
        {
            if(i/n == I || i%n == J)
            {}
            else
            {
                Minor.data[count] = data[i];
                count ++;
            }
        }
        return Minor;
    }

    virtual Matrix reverse()
    {
        if(n!=m || failed())
            return *this;
        float det = determinant();
        if (det == 0)
            return *this;
        Matrix reverse(n,n);
        for(int i=0; i<n; i++)
            for(int j = 0; j<n; j++)
            {
                Matrix Minor = minor(i, j);
                float d = Minor.determinant()/det;
                reverse.set(i,j,d);
                det = -det;
            }
        return reverse;
    }
    virtual Matrix transpose()
    {
        if( failed() )
        {
            Matrix trans;
            return trans;
        }
        Matrix trans(n, m);
        for(int i=0; i<m; i++)
            for(int j=0; j<n; j++)
                trans.set(j,i,my_get(i,j));
        return trans;
    }
    virtual float determinant()
    {
        if(n!=m || failed())
            return 0;
        float Data[n*n];
        for(int i=0; i<n*n; i++)
            Data[i] = data[i];
        float det = determinant(Data, n);
        return det;
    }
    virtual ostream& print(ostream& o)
    {
        for(int i=0; i<m; i++)
        {
            for(int j=0; j<n; j++)
            {
                o << this->my_get(i,j) << '\t';
            }
            o << endl;
        }
        return o;
    }
    virtual istream& read(istream& o)
    {
        o >> this->m >> this->n;
        if(data!=NULL)
            delete data;
        data = new float [m*n];
        for(int i=0; i<m*n; i++)
            o >> data[i];
        return o;
    }
    virtual void set(int i, int j, float data)
    {
        this->data[i*n+j] = data;
        return;
    }
    virtual float my_get(int i, int j)
    {
        return data[i*n+j];
    }
    virtual float get(int i, int j)
    {
        return data[(i-1)*n+(j-1)];
    }
    virtual int getN(){return n;}
    virtual int getM(){return m;}
    virtual bool failed() {return (n<=0 || m<=0 || n!=m || data == NULL);}
};

Matrix* get_init(int n, int m)
{
    Matrix* N = new Matrix(m,n);
    return N;
}

 /*

void foo(Matrix& t)
{
    for(;;)
    t.reverse();
}

int main()
{
    int n = 3;
    Matrix t(n,n);
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            t.set(i,j,1);
    for(int i=0; i<n; i++)
        t.set(i,i,2);
    t.print(cout);
    cout << t.determinant() << endl;
    foo(t);
    return 0;
}*/


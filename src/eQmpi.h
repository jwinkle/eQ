#ifndef EQMPI_H
#define EQMPI_H

#include <mpi.h>
#include <memory>
#include <vector>

namespace eQ {


//================================================================================
//                  MPI TRANSFER UTILITY CLASS
//================================================================================
class mpi
{
public:
	enum class method
	{
		SEND, RECV, ISEND, IRECV,
		BROADCAST,
		WAIT_FOR_ACKNOWLEDGE,
		IDLE
	};

	virtual ~mpi()	=default;
	mpi()			=default;
	//init with node that calling node communicates with (one-to-one) on comm
	void init(MPI_Comm comm, int nodeNumber)
	{
		_comm = comm;
		_nodeNumber = nodeNumber;
		_initialized = true;
		//used by node 0 for indexing vectors of size of diffusion-layer from 0:
		_index = (nodeNumber > 0) ? (nodeNumber - 1) : 0;
		resetTransfer();
	}
	size_t  index(){return _index;}
	bool    waiting(){return _waiting;}

//	mpi &operator>>(std::shared_ptr<std::vector<double>> data);
//	mpi &operator<<(std::shared_ptr<std::vector<double>> data);
//	mpi &operator>>(std::shared_ptr<std::vector<int>> data);
//	mpi &operator<<(std::shared_ptr<std::vector<int>> data);
//	mpi &operator>>(std::vector<double> &data);
//	mpi &operator<<(std::vector<double> &data);
//	mpi &operator>>(std::vector<int> &data);
//	mpi &operator<<(std::vector<int> &data);
//	mpi &operator>>(mpi::method onTask);
//	mpi &operator<<(mpi::method onTask);
protected:
	void setWait() {_waiting = true;}
	void resetWait() {_waiting = false;}
	void setState(mpi::method state) {_mpiState = state;}
	void resetTransfer()
	{
		resetWait();
		mpiRequest.clear();
		setState(mpi::method::IDLE);
	}
	void initTransfer(mpi::method whichMethod)
	{
		resetTransfer();
		setState(whichMethod);
	}

	bool            _initialized=false;
	bool            _waiting=false;
	int             _nodeNumber;
	size_t          _index;
	MPI_Comm        _comm;
	mpi::method       _mpiState=mpi::method::IDLE;

	std::vector<MPI_Request>     mpiRequest;

public:
	//================================================================================
	//                      CONTROL METHOD
	//================================================================================
	//                      RECEIVE, WAIT
		mpi &
		operator>>(mpi::method onTask)
		{
			MPI_Status thisStaus;
			switch(onTask)
			{
			case mpi::method::BROADCAST:
				initTransfer(onTask);
				break;
			case mpi::method::IRECV:
				initTransfer(onTask);
				break;
			case mpi::method::RECV:
				initTransfer(onTask);
				break;
			  case mpi::method::WAIT_FOR_ACKNOWLEDGE:
				if(_waiting)
				{
					for(auto &request : mpiRequest)
					{
						MPI_Wait(&request, &thisStaus);
					}
				}
				resetTransfer();
				break;
			  default:
				resetTransfer();
				break;
			}
			return *this;
		}
	//                      SEND
		mpi &
		operator<<(mpi::method onTask)
		{
			switch(onTask)
			{
				case mpi::method::BROADCAST:
					initTransfer(onTask);
					break;
				case mpi::method::ISEND:
					initTransfer(onTask);
					break;
				case mpi::method::SEND:
					initTransfer(onTask);
					break;
				default:
				  resetTransfer();
					break;
			}
			return *this;
		}
	//================================================================================
	//                  BROADCAST	long
	//================================================================================
	mpi &//				ALL "RECEIVE" BROADCAST DATA (but _nodeNumber actually sends)
	operator>>(long  &data)
	{
		int mpiTag=0;
		MPI_Status thisStaus;

		if(mpi::method::BROADCAST == _mpiState)
		{
			MPI_Bcast(&data, 1, MPI_LONG,
					 _nodeNumber, _comm);
		}
	}
	//================================================================================
	//                      shared_ptr<vector<double>>
	//================================================================================
	//                  RECEIVE
	mpi &
	operator>>(std::shared_ptr<std::vector<double>> data)
	{
		int mpiTag=0;
		MPI_Status thisStaus;

		if(mpi::method::RECV == _mpiState)
		{
			MPI_Recv(data->data(), int(data->size()), MPI_DOUBLE,
					 _nodeNumber, mpiTag, _comm, &thisStaus);
		}
		else if(mpi::method::IRECV == _mpiState)
		{
			MPI_Request request;
			mpiRequest.push_back(request);
			MPI_Irecv(data->data(), int(data->size()), MPI_DOUBLE,
					  _nodeNumber, mpiTag, _comm, &mpiRequest.back());
			setWait();
		}
		return *this;
	}
	//                      SEND
		mpi &
		operator<<(std::shared_ptr<std::vector<double>> data)
		{
			int mpiTag=0;

			if(mpi::method::SEND == _mpiState)
			{
				MPI_Send(data->data(), int(data->size()), MPI_DOUBLE,
						 _nodeNumber, mpiTag, _comm);
			}
			else if(mpi::method::ISEND == _mpiState)
			{
				MPI_Request request;
				mpiRequest.push_back(request);
				MPI_Isend(data->data(), int(data->size()), MPI_DOUBLE,
						  _nodeNumber, mpiTag, _comm, &mpiRequest.back());
				setWait();
			}
			return *this;
	}
	//================================================================================
	//                      shared_ptr<vector<int>>
	//================================================================================
	//                      RECEIVE
		mpi &
		operator>>(std::shared_ptr<std::vector<int>> data)
		{
			int mpiTag=0;
			MPI_Status thisStaus;

			if(mpi::method::RECV == _mpiState)
			{
				MPI_Recv(data->data(), int(data->size()), MPI_INT,
						 _nodeNumber, mpiTag, _comm, &thisStaus);
			}
			else if(mpi::method::IRECV == _mpiState)
			{
				MPI_Request request;
				mpiRequest.push_back(request);
				MPI_Irecv(data->data(), int(data->size()), MPI_INT,
						  _nodeNumber, mpiTag, _comm, &mpiRequest.back());
				setWait();
			}
			return *this;
		}
	//                      SEND
		mpi &
		operator<<(std::shared_ptr<std::vector<int>> data)
		{
			int mpiTag=0;

			if(mpi::method::SEND == _mpiState)
			{
				MPI_Send(data->data(), int(data->size()), MPI_INT,
						 _nodeNumber, mpiTag, _comm);
			}
			else if(mpi::method::ISEND == _mpiState)
			{
				MPI_Request request;
				mpiRequest.push_back(request);
				MPI_Isend(data->data(), int(data->size()), MPI_INT,
						  _nodeNumber, mpiTag, _comm, &mpiRequest.back());
				setWait();
			}
			return *this;
		}
	//================================================================================
	//                      <vector<double> &
	//================================================================================
	//                      RECEIVE
		mpi &
		operator>>(std::vector<double> &data)
		{
			int mpiTag=0;
			MPI_Status thisStaus;

			if(mpi::method::RECV == _mpiState)
			{
				MPI_Recv(data.data(), int(data.size()), MPI_DOUBLE,
						 _nodeNumber, mpiTag, _comm, &thisStaus);
			}
			else if(mpi::method::IRECV == _mpiState)
			{
				MPI_Request request;
				mpiRequest.push_back(request);
				MPI_Irecv(data.data(), int(data.size()), MPI_DOUBLE,
						  _nodeNumber, mpiTag, _comm, &mpiRequest.back());
				setWait();
			}
			return *this;
		}
	//                      SEND
		mpi &
		operator<<(std::vector<double> &data)
		{
			//ship the new HSL grids to the controller node:
			int mpiTag=0;

			if(mpi::method::SEND == _mpiState)
			{
				MPI_Send(data.data(), int(data.size()), MPI_DOUBLE,
						 _nodeNumber, mpiTag, _comm);
			}
			else if(mpi::method::ISEND == _mpiState)
			{
				MPI_Request request;
				mpiRequest.push_back(request);
				MPI_Isend(data.data(), int(data.size()), MPI_DOUBLE,
						  _nodeNumber, mpiTag, _comm, &mpiRequest.back());
				setWait();
			}
			return *this;
		}
	//================================================================================
	//                      <vector<int> &
	//================================================================================
	//                      RECEIVE
		mpi &
		operator>>(std::vector<int> &data)
		{
			int mpiTag=0;
			MPI_Status thisStaus;

			if(mpi::method::RECV == _mpiState)
			{
				MPI_Recv(data.data(), int(data.size()), MPI_INT,
						 _nodeNumber, mpiTag, _comm, &thisStaus);
			}
			else if(mpi::method::IRECV == _mpiState)
			{
				MPI_Request request;
				mpiRequest.push_back(request);
				MPI_Irecv(data.data(), int(data.size()), MPI_INT,
						  _nodeNumber, mpiTag, _comm, &mpiRequest.back());
				setWait();
			}
			return *this;
		}
	//                      SEND
		mpi &
		operator<<(std::vector<int> &data)
		{
			//ship the new HSL grids to the controller node:
			int mpiTag=0;

			if(mpi::method::SEND == _mpiState)
			{
				MPI_Send(data.data(), int(data.size()), MPI_INT,
						 _nodeNumber, mpiTag, _comm);
			}
			else if(mpi::method::ISEND == _mpiState)
			{
				MPI_Request request;
				mpiRequest.push_back(request);
				MPI_Isend(data.data(), int(data.size()), MPI_INT,
						  _nodeNumber, mpiTag, _comm, &mpiRequest.back());
				setWait();
			}
			return *this;
		}

	//================================================================================
	//      END     mpi  definitions
	//================================================================================

};

}
#endif // EQMPI_H

#include "threads.hpp"

paz::LockedCv::LockedCv() : _mx(std::make_shared<std::mutex>()), _cv(std::
    make_shared<std::condition_variable>()), _done(std::make_shared<bool>(
    false)) {}

void paz::LockedCv::release()
{
    std::lock_guard lk(*_mx);
    *_done = true;
    _cv->notify_all();
}

void paz::LockedCv::wait() const
{
    std::unique_lock lk(*_mx);
    _cv->wait(lk, [this](){ return *_done; });
}

paz::Threadpool::Threadpool() : _done(false)
{
    const unsigned int numThreads = std::max(1u, std::thread::
        hardware_concurrency());
    _threads.reserve(numThreads);
    for(std::size_t i = 0; i < numThreads; ++i)
    {
        _threads.emplace_back([this]()
        {
            while(true)
            {
                std::unique_lock lk(_mx);
                _cv.wait(lk, [this](){ return _done || !_tasks.empty(); });
                if(_done)
                {
                    return;
                }
                auto task = _tasks.front();
                _tasks.pop();
                lk.unlock();
                task.first();
                task.second.release();
            }
        }
        );
    }
}

paz::Threadpool::~Threadpool()
{
    {
        std::lock_guard lk(_mx);
        _done = true;
        _cv.notify_all();
    }
    for(auto& n : _threads)
    {
        n.join();
    }
}

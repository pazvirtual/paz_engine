#ifndef PAZ_ENGINE_THREADS_HPP
#define PAZ_ENGINE_THREADS_HPP

#include <thread>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <vector>
#include <functional>
#include <queue>

namespace paz
{
    class LockedCv
    {
        friend class Threadpool;

        std::shared_ptr<std::mutex> _mx;
        std::shared_ptr<std::condition_variable> _cv;
        std::shared_ptr<bool> _done;

        void release();

    public:
        LockedCv();
        void wait() const;
    };

    class Threadpool
    {
        std::mutex _mx;
        std::condition_variable _cv;
        bool _done;
        std::queue<std::pair<std::function<void(void)>, LockedCv>> _tasks;
        std::vector<std::thread> _threads;

    public:
        Threadpool();
        ~Threadpool();
        template<typename F>
        LockedCv pushTask(const F& fun)
        {
            LockedCv lcv;
            std::lock_guard lk(_mx);
            _tasks.emplace(fun, lcv);
            _cv.notify_one();
            return lcv;
        }
    };
}

#endif

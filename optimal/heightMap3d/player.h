#ifndef PLAYER_H
#define PLAYER_H

#include <QObject>
#include <QTimer>

class Player : public QTimer
{
    Q_OBJECT
public:
    explicit Player(QObject *parent = nullptr);

signals:

public slots:
    void timeout();
};

#endif // PLAYER_H

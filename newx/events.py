'''Publish/Subscribe. Preferred way of message sending in newx

    We simply use the base supplied in wx, with some name changes for my
    liking. Mainly: 
        - the topic is a string with '.' separators (like the logging system).
        - the order is consistent for publish and subscribe (topic first)

'''

from wx.lib.pubsub import Publisher
publisher = Publisher()

def publish(path, data=None):
    assert isinstance(path, str)
    topic = tuple(path.split('.'))
    # print topic, data
    publisher.sendMessage(topic, data)

def subscribe(path, listener):
    assert isinstance(path, str)
    topic = tuple(path.split('.'))
    # print listener, topic
    publisher.subscribe(listener, topic)

def unsubscribe(path, listener):
    assert isinstance(path, str)
    topic = tuple(path.split('.'))
    publisher.unsubscribe(listener, topic)

# TODO For this to be useful, it should really be in the base of everything
# class Subscriber(object):
    #
    # '''Provides callback for anyone wanting to subscribe

    # '''
    # def __init__(self):
        # self.make_subscriptions()

    # def make_subscriptions(self):
        # pass


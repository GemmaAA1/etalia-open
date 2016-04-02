define(['lib/backbone'], function(Backbone) {

    var _sync = Backbone.sync;
    Backbone.sync = function(method, model, options){

        // Add trailing slash to backbone model views
        var _url = _.isFunction(model.url) ?  model.url() : model.url;
        _url += _url.charAt(_url.length - 1) == '/' ? '' : '/';

        options = _.extend(options, {
            url: _url
        });

        return _sync(method, model, options);
    };

    return Backbone;
});

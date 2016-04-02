define([
    'app/config',
    'backbone',
    'model/thread/thread'
], function (Config, Backbone, ThreadModel) {

    return Backbone.Collection.extend({
        url: Config.api_root + '/thread/threads',
        model: ThreadModel,

        parse: function(response) {
            return response.results;
        }
    });

});

define([
    'app/config',
    'backbone',
    'model/paper'
], function (Config, Backbone, PaperModel) {

    return Backbone.Collection.extend({
        url: Config.api_root + '/papers',
        model: PaperModel,

        parse: function(response) {
            return response.results;
        }
    });

});

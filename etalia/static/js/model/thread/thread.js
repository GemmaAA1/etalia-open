define([
    'app/config',
    'backbone',
    'model/user/user',
    'model/paper',
    'model/user/library',
    'backbone-relational'
], function (Config, Backbone, UserModel, PaperModel, LibraryModel) {

    var defaults = {
        user: null,
        paper: null,
        type: 2,
        privacy: 1,
        title: null,
        content: null,
        state: null,
        members: [],
        posts: [],
        created: null,
        modified: null
    };

    var schema = {
        title: {type: 'Text', validators: ['required']},
        paper: {type: 'Select', options: function(callback) {

            function parseLibrary(library, callback) {
                var library_papers = library.get('papers'),
                    papers = new Backbone.Collection(
                        library_papers.pluck('paper'),
                        {model: PaperModel}
                    );

                console.log(papers);

                callback(papers);
            }

            var user = UserModel.getCurrent(),
                library = LibraryModel.find(user.get('id'));

            // If library found
            if (library) {
                parseLibrary(library, callback);
                return;
            }

            // Fetch the user library
            library = new LibraryModel({pk: user.get('id')});
            library.fetch().then(function() {
                parseLibrary(library, callback);
            }, function(jqXHR, textStatus, errorThrown) {
                throw textStatus;
            });

        }, validators: ['required']}
    };

    var ThreadModel = Backbone.RelationalModel.extend({
        urlRoot: Config.api_root + '/thread/threads',

        defaults: defaults,

        schema: schema,

        relations: [
            {
                type: Backbone.HasOne,
                key: 'user',
                relatedModel: UserModel,
                includeInJSON: 'link'
            },
            {
                type: Backbone.HasOne,
                key: 'paper',
                relatedModel: PaperModel,
                includeInJSON: 'link'
            }
        ]
    });

    ThreadModel.createNew = function() {
        return new ThreadModel({
            'user': UserModel.getCurrent()
        });
    };

    return ThreadModel;
});
